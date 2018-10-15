import core.stdc.stdio;
//import core.stdc.stdlib;
import std.algorithm.comparison;
import std.file;
import std.experimental.allocator;
import std.stdio: writeln;

/*import <string>
  import <vector>
  import <map>
*/
import read;
import options;
import util;

enum KMER_LEN = 5;

// TODO: since this module uses a lot of string concat
// we might want to use our own string_builder
// which would preallocate

T[] allocArray (T, size_t line = __LINE__, string file = __FILE__) (size_t count)
{
  import std.conv : to;
  import core.stdc.stdlib : calloc;

  static assert(!is(typeof(T.init[0])),
                "T should not be an array -- called from" ~ file ~ ":" ~ line.to!string);

  T* memory = cast(T*) calloc(count, T.sizeof);
  auto result = memory[0 .. count];
  result[] = T.init;
  return result;
}

void free(void* mem/*, size_t line = __LINE__, string file = __FILE__*/)
{
  import core.stdc.stdlib;
  core.stdc.stdlib.free(mem);
}


class Stats
{
public:
  // this @guessedCycles parameter should be calculated using the first several records
  this(Options* opt, bool isRead2 = false, size_t guessedCycles = 0, size_t bufferMargin = 1024)
  {
    mOptions = opt;
    mIsRead2 = isRead2;
    mReads = 0;

    mEvaluatedSeqLen = mOptions.seqLen1;
    if (mIsRead2)
      mEvaluatedSeqLen = mOptions.seqLen2;

    if (guessedCycles == 0)
      {
        guessedCycles = mEvaluatedSeqLen;
      }

    mCycles = guessedCycles;
    mBases = 0;
    mQ20Total = 0;
    mQ30Total = 0;
    summarized = false;
    mKmerMin = 0;
    mKmerMax = 0;

    // extend the buffer to make sure it's long enough
    mBufLen = cast(size_t)(guessedCycles + bufferMargin);

    //TODO: Repetitive code below, refactoring needed
    foreach (i ; 0 .. 8 )
      {
        mQ20Bases[i] = mQ30Bases[i] = mBaseContents[i] = 0;
        // mCycleQ30Bases[mBufLen][8]
        // Allocate an array of mBufLen number of longs and initialize all of them to 0
        mCycleQ30Bases[i] = allocArray!long(mBufLen);
        mCycleQ20Bases[i] = allocArray!long(mBufLen);
        mCycleBaseContents[i] = allocArray!long(mBufLen);
        mCycleBaseQual[i] = allocArray!long(mBufLen);
      }

    mCycleTotalBase = allocArray!long(mBufLen);
    mCycleTotalQual= allocArray!long(mBufLen);
    mKmerBufLen = 2 << (KMER_LEN * 2);
    mKmer = allocArray!long(mKmerBufLen);
    initOverRepSeq();
  }

  ~this()
  {
    foreach (i ; 0 .. 8 )
      {
        free(mCycleQ30Bases[i].ptr);
        free(mCycleQ20Bases[i].ptr);
        free(mCycleBaseContents[i].ptr);
        free(mCycleBaseQual[i].ptr);
      }

    free(mCycleTotalBase.ptr);
    free(mCycleTotalQual.ptr);

    // delete memory of curves

    // Note: We may call garbage collector actively here if needed.
    mQualityCurves.destroy;
    mContentCurves.destroy;
    free(&mKmer[0]); //The compiler does not complain when we use @safe;
    deleteOverRepSeqDist();
  }

  size_t getCycles()
  {
    if (!summarized)
      summarize();
    return mCycles;
  }

  long getReads()
  {
    if (!summarized)
      summarize();
    return mReads;
  }

  long getBases()
  {
    if (!summarized)
      summarize();
    return mBases;
  }

  long getQ20()
  {
    if (!summarized)
      summarize();
    return mQ20Total;
  }

  long getQ30()
  {
    if (!summarized)
      summarize();
    return mQ30Total;
  }

  long getGCNumber()
  {
    if (!summarized)
      summarize();
    return mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07];
  }
  // by default the qualified qual score is Q20 ('5')
  void statRead(Read r)
  {
    size_t len = r.length;

    if (mBufLen < len)
      {
        extendBuffer(max(len + 100, cast(size_t)(len * 1.5)));
      }
    const (char)[] seqstr = r.mSeq.mStr;
    const (char)[] qualstr = r.mQuality;

    int kmer = 0;
    bool needFullCompute = true;
    for (int i = 0; i < len; i++)
      {
        char base = seqstr[i];
        char qual = qualstr[i];
        // get last 3 bits
        char b = base & 0x07;

        const char q20 = '5';
        const char q30 = '?';

        if (qual >= q30)
          {
            mCycleQ30Bases[b][i]++;
            mCycleQ20Bases[b][i]++;
          }
        else if (qual >= q20)
          {
            mCycleQ20Bases[b][i]++;
          }

        mCycleBaseContents[b][i]++;
        mCycleBaseQual[b][i] += (qual - 33);

        mCycleTotalBase[i]++;
        mCycleTotalQual[i] += (qual - 33);

        if (base == 'N')
          {
            needFullCompute = true;
            continue;
          }

        // 5 bases required for kmer computing
        if (i < 4)
          continue;

        // calc 5 KMER
        // 0x3FC == 0011 1111 1100
        if (!needFullCompute)
          {
            int val = base2val(base);
            if (val < 0)
              {
                needFullCompute = true;
                continue;
              }
            else
              {
                kmer = ((kmer << 2) & 0x3FC) | val;
                mKmer[kmer]++;
              }
          }
        else
          {
            bool valid = true;
            kmer = 0;
            for (int k = 0; k < 5; k++)
              {
                int val = base2val(seqstr[i - 4 + k]);
                if (val < 0)
                  {
                    valid = false;
                    break;
                  }
                kmer = ((kmer << 2) & 0x3FC) | val;
              }
            if (!valid)
              {
                needFullCompute = true;
                continue;
              }
            else
              {
                mKmer[kmer]++;
                needFullCompute = false;
              }
          }

      }

    // do overrepresentation analysis for 1 of every 100 reads
    if (mOptions.overRepAnalysis.enabled)
      {
        if (mReads % mOptions.overRepAnalysis.sampling == 0)
          {
            const int[5] steps = [10, 20, 40, 100, min(150, mEvaluatedSeqLen - 2)];
            for (int s = 0; s < 5; s++)
              {
                int step = steps[s];
                for (int i = 0; i < len - step; i++)
                  {
                    string seq = r.mSeq.mStr[i .. i + step];
                    if (auto elem = seq in mOverRepSeq)
                      {
                        (*elem)++;
                        for (int p = i; p < seq.length + i && p < mEvaluatedSeqLen;
                             p++)
                          {
                            mOverRepSeqDist[seq][p]++;
                          }
                        i += step;
                      }
                  }
              }
          }
      }

    mReads++;
  }

  private Stats merge(Stats[] list)
  {
    if (list.length == 0)
      return null;

    //get the most long cycles
    size_t cycles = 0;
    for (size_t t = 0; t < list.length; t++)
      {
        list[t].summarize();
        cycles = max(cycles, list[t].getCycles());
      }

    Stats s = new Stats(list[0].mOptions, list[0].mIsRead2, cycles, 0);

    // init overrepresented seq maps
    long[string] iter;

    for (size_t t = 0; t < list.length; t++)
      {
            size_t curCycles = list[t].getCycles();
        // merge read number
        s.mReads += list[t].mReads;

        // merge per cycle counting for different bases
        for (int i = 0; i < 8; i++)
          {
            for (int j = 0; j < cycles && j < curCycles; j++)
              {
                s.mCycleQ30Bases[i][j] += list[t].mCycleQ30Bases[i][j];
                s.mCycleQ20Bases[i][j] += list[t].mCycleQ20Bases[i][j];
                s.mCycleBaseContents[i][j] += list[t].mCycleBaseContents[i][j];
                s.mCycleBaseQual[i][j] += list[t].mCycleBaseQual[i][j];
              }
          }

        // merge per cycle counting for all bases
        for (int j = 0; j < cycles && j < curCycles; j++)
          {
            s.mCycleTotalBase[j] += list[t].mCycleTotalBase[j];
            s.mCycleTotalQual[j] += list[t].mCycleTotalQual[j];
          }

        // merge kMer
        for (int i = 0; i < s.mKmerBufLen; i++)
          {
            s.mKmer[i] += list[t].mKmer[i];
          }

        // merge over rep seq
        foreach (string seq,  long c ; s.mOverRepSeq)
          {
            s.mOverRepSeq[seq] += list[t].mOverRepSeq[seq];
            if (s.mIsRead2 != list[t].mIsRead2
                || list[t].mOverRepSeqDist[seq] == null)
                {
                writeln(
                  t, seq, ":", (s.mIsRead2 ? 2 : 1), ",",
                  (list[t].mIsRead2  ? 2 : 1)
                );
              }
          for (int i = 0; i < s.mEvaluatedSeqLen; i++)
            {
              s.mOverRepSeqDist[seq][i] += list[t].mOverRepSeqDist[seq][i];
            }
        }
    }

  s.summarize();

  return s;
}

void print()
{
  if (!summarized)
    {
      summarize();
    }
    writeln("total reads: ", mReads);
    writeln( "total bases: ", mBases);
    writeln("Q20 bases: " , mQ20Total , "(", (mQ20Total * 100.0) / mBases , "%)");
    writeln("Q30 bases: ", mQ30Total, "(", (mQ30Total * 100.0) / mBases, "%)");
  }

  void summarize(bool forced = false)
  {
    if (summarized && !forced)
      return;

    // first get the cycle and count total bases
    for (size_t c = 0; c < mBufLen; c++)
      {
        mBases += mCycleTotalBase[c];
        if (mCycleTotalBase[c] == 0)
          {
            mCycles = c;
            break;
          }
      }
    if (mCycleTotalBase[mBufLen - 1] > 0)
      mCycles = mBufLen;

    // TODO: Can we shorten the nested for block below?
    // Q20, Q30, base content
    for (int i = 0; i < 8; i++)
      {
        for (int c = 0; c < mCycles; c++)
          {
            mQ20Bases[i] += mCycleQ20Bases[i][c];
            mQ30Bases[i] += mCycleQ30Bases[i][c];
            mBaseContents[i] += mCycleBaseContents[i][c];
          }
        mQ20Total += mQ20Bases[i];
        mQ30Total += mQ30Bases[i];
      }

    // quality curve for mean qual
    double[] meanQualCurve = new double[] (mCycles);
    meanQualCurve[] = 0;
    // TODO: Check and make sure the array operation works as expected.
    meanQualCurve[] =
      (cast(double[])(mCycleTotalQual))[] /
      (cast(double[])(mCycleTotalBase))[];

    mQualityCurves["mean"] = meanQualCurve;

    // quality curves and base content curves for different nucleotides
    static immutable char[5] alphabets = "ATCGN";// can avoid array allocation with this
    for (int i = 0; i < 5; i++)
      {
        char base = alphabets[i];
        // get last 3 bits
        char b = base & 0x07;
        double[] qualCurve = new double[mCycles];
        qualCurve[] = 0;
        double[] contentCurve = new double[mCycles];
        contentCurve[] = 0;
        // TODO: Shorten the for loop
        for (int c = 0; c < mCycles; c++)
          {
            if (mCycleBaseContents[b][c] == 0)
              qualCurve[c] = meanQualCurve[c];
            else
              qualCurve[c] = cast(double) mCycleBaseQual[b][c] / cast(double) mCycleBaseContents[b][c];
            contentCurve[c] = cast(double) mCycleBaseContents[b][c] / cast(double) mCycleTotalBase[c];
          }
        static char[1] base_string;
        base_string[0] = base;
        mQualityCurves[cast(string)base_string] = qualCurve;
        mContentCurves[cast(string)base_string] = contentCurve;
      }

    // GC content curve
    double[] gcContentCurve = new double[mCycles];
    double[] gcContentCurve_tmp = new double[mCycles];
    gcContentCurve[] = 0;
    char gBase = 'G' & 0x07;
    char cBase = 'C' & 0x07;

    gcContentCurve_tmp[] =
        (cast(double[])mCycleBaseContents[gBase])[] +
        (cast(double[])mCycleBaseContents[cBase])[];

    gcContentCurve[] = gcContentCurve_tmp[] /
       (cast(double[])mCycleTotalBase)[];

    mContentCurves["GC"] = gcContentCurve;

    mKmerMin = mKmer[0];
    mKmerMax = mKmer[0];
    for (int i = 0; i < mKmerBufLen; i++)
      {
        if (mKmer[i] > mKmerMax)
          mKmerMax = mKmer[i];
        if (mKmer[i] < mKmerMin)
          mKmerMin = mKmer[i];
      }

    summarized = true;
  }
  // a port of JSON report
  void reportJson(ref char[] outa, string padding)
  {
    outa ~= "{\n";
    import std.conv : to;
    outa ~= padding ~ "\t\"total_reads\": " ~ to!string(mReads) ~ ",\n";
    outa ~= padding ~ "\t\"total_bases\": " ~ to!string(mBases) ~ ",\n";
    outa ~= padding ~ "\t\"q20_bases\": " ~ to!string(mQ20Total) ~ ",\n";
    outa ~= padding ~ "\t\"q30_bases\": " ~ to!string(mQ30Total) ~ ",\n";
    outa ~= padding ~ "\t\"total_cycles\": " ~ to!string(mCycles) ~ ",\n";

    // quality curves
    static immutable string[5] qualNames = ["A", "T", "C", "G", "mean"];
    outa ~= padding ~ "\t" ~ "\"quality_curves\": {" ~ "\n";
    for (int i = 0; i < 5; i++)
      {
        string name = qualNames[i];
        double[] curve = mQualityCurves[name];
        outa ~= padding ~ "\t\t" ~ "\"" ~ name ~ "\":[";
        for (int c = 0; c < mCycles; c++)
          {
            outa ~= to!string(curve[c]);
            // not the end
            if (c != mCycles - 1)
              outa ~= ",";
          }
        outa ~= "]";
        // not the end;
        if (i != 5 - 1)
          outa ~= ",";
        outa ~= "\n";
      }
    outa ~= padding ~ "\t" ~ "}," ~ "\n";

    // content curves
    static immutable string[6] contentNames = ["A", "T", "C", "G", "N", "GC"];
    outa ~= padding ~ "\t" ~ "\"content_curves\": {" ~ "\n";
    for (int i = 0; i < 6; i++)
      {
        string name = contentNames[i];
        double[] curve = mContentCurves[name];
        outa ~= padding ~ "\t\t" ~ "\"" ~ name ~ "\":[";
        for (int c = 0; c < mCycles; c++)
          {
            outa ~= to!string(curve[c]);
            // not the end
            if (c != mCycles - 1)
              outa ~= ",";
          }
        outa ~= "]";
        // not the end;
        if (i != 6 - 1)
          outa ~= ",";
        outa ~= "\n";
      }
    outa ~= padding ~ "\t" ~ "}," ~ "\n";

    // KMER counting
    outa ~= padding ~ "\t" ~ "\"kmer_count\": {" ~ "\n";
    for (int i = 0; i < 64; i++)
      {
        string first = kmer3(i);
        for (int j = 0; j < 16; j++)
          {
            int target = (i << 4) + j;
            long count = mKmer[target];
            string last = kmer2(j);
            outa ~= padding ~ "\t\t\"" ~ first ~ last ~ "\":" ~ to!string(count);
            if (j != 16 - 1)
              outa ~= ",";
          }
        if (i != 64 - 1)
          outa ~= "," ~ "\n";
        else
          outa ~= "\n";
      }
    outa ~= padding ~ "\t" ~ "}," ~ "\n";

    // over represented seqs
    long [string] iter;
    bool first = true;
    outa ~= padding ~ "\t" ~ "\"overrepresented_sequences\": {" ~ "\n";
    foreach (seq, count; mOverRepSeq)
      {
        //string seq = iter.first;
        //long count = iter.second;
        if (!overRepPassed(seq, count))
          continue;
        if (!first)
          {
            outa ~= "," ~ "\n";
          }
        else
          first = false;
        outa ~= padding ~ "\t\t\"" ~ seq ~ "\":" ~ to!string(count);
      }
    outa ~= padding ~ "\t" ~ "}" ~ "\n";

    outa ~= padding ~ "}," ~ "\n";
  }

  // a port of HTML report
  void reportHtml(ref char[] outa, string filteringType, string readName)
  {
  // There can be a performance problem because this allocates a lot.
    reportHtmlQuality(outa, filteringType, readName);
    reportHtmlContents(outa, filteringType, readName);
    reportHtmlKMER(outa, filteringType, readName);
    if (mOptions.overRepAnalysis.enabled)
      {
        reportHtmlORA(outa, filteringType, readName);
      }
  }
  void reportHtmlQuality(ref char[] outa, string filteringType, string readName)
  {

    // quality
    string subsection = filteringType ~ ": " ~ readName ~ ": quality";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    outa ~= "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" ~ divName ~ "')>"
      ~ subsection ~ "</a></div>\n";
    outa ~= "<div id='" ~ divName ~ "'>\n";
    outa ~= "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    outa ~= "<div class='figure' id='plot_" ~ divName ~ "'></div>\n";
    outa ~= "</div>\n";

    static immutable string[5] alphabets = ["A", "T", "C", "G", "mean"];
    static immutable string[5] colors = [
      "rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)",
      "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"
    ];
    outa ~= "\n<script type=\"text/javascript\">" ~ "\n";
    string json_str = "var data=[";

    long[] x = new long[mCycles];
    int total = 0;
    if (!isLongRead())
      {
        for (int i = 0; i < mCycles; i++)
          {
            x[total] = i + 1;
            total++;
          }
      }
    else
      {
        const int fullSampling = 40;
        for (int i = 0; i < fullSampling && i < mCycles; i++)
          {
            x[total] = i + 1;
            total++;
          }
        // down sampling if it's too long
        if (mCycles > fullSampling)
          {
            double pos = fullSampling;
            while (true)
              {
                pos *= 1.05;
                if (pos >= mCycles)
                  break;
                x[total] = cast (int) pos;
                total++;
              }
            // make sure lsat one is contained
            if (x[total - 1] != mCycles)
              {
                x[total] = mCycles;
                total++;
              }
          }
      }
    // four bases
    for (int b = 0; b < 5; b++)
      {
        string base = alphabets[b];
        json_str ~= "{";
        json_str ~= "x:[" ~ list2string(&x[0], total) ~ "],";
        json_str ~= "y:[" ~ list2string(&mQualityCurves[base][0], total, &x[0]) ~ "],";
        json_str ~= "name: '" ~ base ~ "',";
        json_str ~= "mode:'lines',";
        json_str ~= "line:{color:'" ~ colors[b] ~ "', width:1}\n";
        json_str ~= "},";
      }
    json_str ~= "];\n";
    json_str ~= "var layout={title:'" ~ title ~ "', xaxis:{title:'position'";
    // use log plot if it's too long
    if (isLongRead())
      {
        json_str ~= ",type:'log'";
      }
    json_str ~= "}, yaxis:{title:'quality'}};\n";
    json_str ~= "Plotly.newPlot('plot_" ~ divName ~ "', data, layout);\n";

    outa ~= json_str;
    outa ~= "</script>" ~ "\n";

    x.destroy();
  }
  void reportHtmlContents(ref char[] outa, string filteringType, string readName)
  {

    // content
        string subsection = filteringType ~ ": " ~ readName ~ ": base contents";
        string divName = replace(subsection, " ", "_");
        divName = replace(divName, ":", "_");
        string title = "";

        outa ~= "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" ~ divName ~ "')>"
            ~ subsection ~ "</a></div>\n";
        outa ~= "<div id='" ~ divName ~ "'>\n";
        outa ~= "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
        outa ~= "<div class='figure' id='plot_" ~ divName ~ "'></div>\n";
        outa ~= "</div>\n";

        static immutable string[6] alphabets = ["A", "T", "C", "G", "N", "GC"];
        static immutable string[6] colors = [
          "rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)",
          "rgba(0,0,255,1.0)", "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"
        ];
        outa ~= "\n<script type=\"text/javascript\">" ~ "\n";
        string json_str = "var data=[";

        long[] x = allocArray!long(mCycles);
        int total = 0;
        if (!isLongRead())
    {
      for (int i = 0; i < mCycles; i++)
        {
          x[total] = i + 1;
          total++;
        }
    }
  else
    {
      const int fullSampling = 40;
      for (int i = 0; i < fullSampling && i < mCycles; i++)
        {
          x[total] = i + 1;
          total++;
        }
      // down sampling if it's too long
      if (mCycles > fullSampling)
        {
          double pos = fullSampling;
          while (true)
            {
              pos *= 1.05;
              if (pos >= mCycles)
                break;
              x[total] = cast(int) pos;
              total++;
            }
          // make sure lsat one is contained
          if (x[total - 1] != mCycles)
            {
              x[total] = mCycles;
              total++;
            }
        }
    }
  // four bases
  for (int b = 0; b < 6; b++)
    {
      string base = alphabets[b];
      long count = 0;
      if (base.length == 1)
        {
          char _b = base[0] & 0x07;
          static assert(0x07 == 0b111);
          count = mBaseContents[_b];
        }
      else
        {
          count = mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07];
        }
        import std.conv : to;
      string percentage = to!string(cast (double) count * 100.0 / mBases);
      if (percentage.length > 5)
        percentage = percentage[0 .. 5];

      string name = base ~ "(" ~ percentage ~ "%)";

      json_str ~= "{";
      json_str ~= "x:[" ~ list2string(&x[0], total) ~ "],";
      json_str ~= "y:[" ~ list2string(&mContentCurves[base][0], total, &x[0]) ~ "],";
      json_str ~= "name: '" ~ name ~ "',";
      json_str ~= "mode:'lines',";
      json_str ~= "line:{color:'" ~ colors[b] ~ "', width:1}\n";
      json_str ~= "},";
    }
  json_str ~= "];\n";
  json_str ~= "var layout={title:'" ~ title ~ "', xaxis:{title:'position'";
  // use log plot if it's too long
  if (isLongRead())
    {
      json_str ~= ",type:'log'";
    }
  json_str ~= "}, yaxis:{title:'base content ratios'}};\n";
  json_str ~= "Plotly.newPlot('plot_" ~ divName ~ "', data, layout);\n";

  outa ~= json_str;
  outa ~= "</script>" ~ "\n";

  free (x.ptr);
  }
  void reportHtmlKMER(ref char[] outa, string filteringType, string readName)
  {

    // KMER
    string subsection = filteringType ~ ": " ~ readName ~ ": KMER counting";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    outa ~= "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" ~ divName ~ "')>"
      ~ subsection ~ "</a></div>\n";
    outa ~= "<div  id='" ~ divName ~ "'>\n";
    outa ~= "<div class='sub_section_tips'>Darker background means larger counts. The count will be shown on mouse over.</div>\n";
    outa ~= "<table class='kmer_table' style='width:680px;'>\n";
    outa ~= "<tr>";
    outa ~= "<td></td>";
    // the heading row
    for (int h = 0; h < 16; h++)
      outa ~= "<td style='color:#333333'>" ~ kmer2(h) ~ "</td>";
    outa ~= "</tr>\n";
    // content
    for (int i = 0; i < 64; i++)
      {
        outa ~= "<tr>";

        outa ~= "<td style='color:#333333'>" ~ kmer3(i) ~ "</td>";
        for (int j = 0; j < 16; j++)
          {
            outa ~= makeKmerTD(i, j);
          }
        outa ~= "</tr>\n";
      }
    outa ~= "</table>\n";
    outa ~= "</div>\n";
  }
  void reportHtmlORA(ref char[] outa, string filteringType, string readName)
  {
    // over represented seqs
    double dBases = mBases;
    long[string] iter;
    int displayed = 0;

    // KMER
    string subsection = filteringType ~ ": " ~ readName ~ ": overrepresented sequences";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    outa ~= "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" ~ divName ~ "')>"
      ~ subsection ~ "</a></div>\n";
    outa ~= "<div  id='" ~ divName ~ "'>\n";
    import std.conv : to;
    outa ~= "<div class='sub_section_tips'>Sampling rate: 1 / " ~
      to!string(mOptions.overRepAnalysis.sampling) ~ "</div>\n";
    outa ~= "<table class='summary_table'>\n";
    outa ~= "<tr style='font-weight:bold;'><td>overrepresented sequence</td><td>count (% of bases)</td><td>distribution: cycle 1 ~ cycle " ~ to!string(mEvaluatedSeqLen) ~ "</td></tr>" ~ "\n";
    int found = 0;
    foreach (seq, count; mOverRepSeq)
      {
        if (!overRepPassed(seq, count))
          continue;
        found++;
        double percent = (100.0 * count * seq.length * mOptions.overRepAnalysis.sampling) / dBases;
        outa ~= "<tr>";
        outa ~= "<td width='400' style='word-break:break-all;font-size:8px;'>" ~ seq ~ "</td>";
        outa ~= "<td width='200'>" ~ to!string(count) ~ " (" ~ to!string(percent) ~ "%)</td>";
        outa ~= "<td width='250'><canvas id='" ~ divName ~ "_" ~ seq ~ "' width='240' height='20'></td>";
        outa ~= "</tr>" ~ "\n";
      }
    if (found == 0)
      outa ~= "<tr><td style='text-align:center' colspan='3'>not found</td></tr>" ~ "\n";
    outa ~= "</table>\n";
    outa ~= "</div>\n";

    // output the JS
    outa ~= "<script language='javascript'>" ~ "\n";
    outa ~= "var seqlen = " ~ to!string(mEvaluatedSeqLen) ~ ";" ~ "\n";
    outa ~= "var orp_dist = {" ~ "\n";
    bool first = true;
    foreach (seq, count; mOverRepSeq)
      {
        if (!overRepPassed(seq, count))
          continue;

        if (!first)
          {
            outa ~= "," ~ "\n";
          }
        else
          first = false;
        outa ~= "\t\"" ~ divName ~ "_" ~ seq ~ "\":[";
        for (int i = 0; i < mEvaluatedSeqLen; i++)
          {
            if (i != 0)
              outa ~= ",";
            outa ~= to!string(mOverRepSeqDist[seq][i]);
          }
        outa ~= "]";
      }
    outa ~= "\n};" ~ "\n";

    outa ~= "for (seq in orp_dist) {" ~ "\n";
    outa ~= "    var cvs = document.getElementById(seq);" ~ "\n";
    outa ~= "    var ctx = cvs.getContext('2d'); " ~ "\n";
    outa ~= "    var data = orp_dist[seq];" ~ "\n";
    outa ~= "    var w = 240;" ~ "\n";
    outa ~= "    var h = 20;" ~ "\n";
    outa ~= "    ctx.fillStyle='#cccccc';" ~ "\n";
    outa ~= "    ctx.fillRect(0, 0, w, h);" ~ "\n";
    outa ~= "    ctx.fillStyle='#0000FF';" ~ "\n";
    outa ~= "    var maxVal = 0;" ~ "\n";
    outa ~= "    for(d=0; d<seqlen; d++) {" ~ "\n";
    outa ~= "        if(data[d]>maxVal) maxVal = data[d];" ~ "\n";
    outa ~= "    }" ~ "\n";
    outa ~= "    var step = (seqlen-1) /  (w-1);" ~ "\n";
    outa ~= "    for(x=0; x<w; x++){" ~ "\n";
    outa ~= "        var target = step * x;" ~ "\n";
    outa ~= "        var val = data[Math.floor(target)];" ~ "\n";
    outa ~= "        var y = Math.floor((val / maxVal) * h);" ~ "\n";
    outa ~= "        ctx.fillRect(x,h-1, 1, -y);" ~ "\n";
    outa ~= "    }" ~ "\n";
    outa ~= "}" ~ "\n";
    outa ~= "</script>" ~ "\n";
  }
  bool isLongRead()
  {
    return mCycles > 300;
  }

  void initOverRepSeq()
  {
    long[string] overRepSeq;
    if (mIsRead2)
      overRepSeq = mOptions.overRepSeqs2;
    else
      overRepSeq = mOptions.overRepSeqs1;

    long[string] iter;
    foreach (seq; mOverRepSeq.keys)
      {
        mOverRepSeq[seq] = 0;
        long[] distBuf = new long[mEvaluatedSeqLen];
        distBuf[] = 0;
        mOverRepSeqDist[seq] = distBuf;
      }
  }
import std.conv : to;
public:
  // TODO: what to do with `static` here?
  static string list2string(double* list, int size)
  {
    import std.conv: to;
    string ss;

    if (size > 0)
        ss ~= to!string(list[0]);

    for (int i = 1; i < size; i++)
      {
        ss ~= ",";
        ss ~= to!string(list[i]);
      }

    return ss;
  }

  static string list2string(double* list, int size, long* coords)
  {
    string ss;
    for (int i = 0; i < size; i++)
      {
        // coords is 1,2,3,...
        long start = 0;
        if (i > 0)
          start = coords[i - 1];
        long end = coords[i];

        double total = 0.0;
        for (long k = start; k < end; k++)
          total += list[k];

        // get average
        if (end == start)
          ss ~= "0.0";
        else
          ss ~= to!string(total / (end - start));
        // ss ~= list[coords[i]-1];
        if (i < size - 1)
          ss ~= ",";
      }
    return ss;
  }

  static string list2string(long* list, int size)
  {
    import std.conv: to;
    string ss;

    if (size > 0)
        ss ~= to!string(list[0]);

    for (int i = 1; i < size; i++)
      {
        ss ~= ",";
        ss ~= to!string(list[i]);
      }

    return ss;
  };

  static int base2val(char base)
  {
    switch (base)
      {
      case 'A':
        return 0;
      case 'T':
        return 1;
      case 'C':
        return 2;
      case 'G':
        return 3;
      default:
        return -1;
      }
  }

private:
  void extendBuffer(size_t newBufLen)
  {
    import core.stdc.string: memcpy, memset;
    if (newBufLen <= mBufLen)
      return;

    long[] newBuf;

    for (int i = 0; i < 8; i++)
      {
               mixin(() {
                 string result;

                 foreach(name; ["mCycleQ30Bases", "mCycleQ20Bases", "mCycleBaseContents", "mCycleBaseQual"])
                   {
                     result ~= q{
                       newBuf = allocArray!long(newBufLen);
                       newBuf[0 .. mBufLen] = } ~ name ~ q{[i][0 .. mBufLen];
                       memcpy(newBuf.ptr, } ~ name ~ q{[i].ptr, long.sizeof * mBufLen);
                       free( } ~ name ~ q{[i].ptr);
                     } ~ name ~ q{[i] = newBuf; }
                     ;
                   }
                 return result;
               } ());

      }

    newBuf = new long[newBufLen];
    memset(&newBuf[0], 0, long.sizeof * newBufLen);
    // &newBuf[0] is better than newBuf.ptr because the compiler can know the
    // length if there is a length
    memcpy(&newBuf[0], &mCycleTotalBase[0], long.sizeof * mBufLen);
    mCycleTotalBase.destroy();
    mCycleTotalBase = newBuf;

    newBuf = new long[newBufLen];
    memset(&newBuf[0], 0, long.sizeof * newBufLen);
    memcpy(&newBuf[0], &mCycleTotalQual[0], long.sizeof * mBufLen);
    mCycleTotalQual.destroy();
    mCycleTotalQual = newBuf;

    mBufLen = newBufLen;
  }

  string makeKmerTD(int i, int j)
  {
    int target = (i << 4) + j;
    long val = mKmer[target];
    // 3bp + 2bp = 5bp
    string first = kmer3(i);
    string last = kmer2(j);
    string kmer = first ~ last;
    double meanBases = cast(double)(mBases + 1) / mKmerBufLen;
    double prop = val / meanBases;
    double frac = 0.5;
    if (prop > 2.0)
      frac = (prop - 2.0) / 20.0 + 0.5;
    else if (prop < 0.5)
      frac = prop;

    frac = max(0.01, min(1.0, frac));
    int r = cast(int) ((1.0 - frac) * 255);
    int g = r;
    int b = r;
    char[] ss;
    ss.length = 255;

    int slen = sprintf(&ss[0], "<td style='background:#%.2X%.2X%.2X 'title='%s: %d\n %d times as mean value'> </td>",
                                             r,  g,  b,           &kmer[0],val, prop);

    return cast(string) ss[0 .. slen];
  }

  string kmer3(int val)
  {
    const char[4] bases = ['A', 'T', 'C', 'G'];
    char[] ret;
    ret.length = 3;
    ret[0] = bases[(val & 0x30) >> 4];
    ret[1] = bases[(val & 0x0C) >> 2];
    ret[2] = bases[(val & 0x03)];
    return cast(string) ret;
  }

  string kmer2(int val)
  {
    const char[4] bases = ['A', 'T', 'C', 'G'];
    char[] ret;
    ret.length = 2;
    ret[0] = bases[(val & 0x0C) >> 2];
    ret[1] = bases[(val & 0x03)];
    return cast(string) ret;
  }

  void deleteOverRepSeqDist()
  {
    // TODO: Deal with loop and delete vector elements
    long[string] iter;
    foreach (seq, _; mOverRepSeq)
      {
        mOverRepSeqDist[seq].destroy();
        mOverRepSeqDist[seq] = null;
      }
  }

  bool overRepPassed(ref string seq, long count)
  {
    int s = mOptions.overRepAnalysis.sampling;
    switch (seq.length)
      {
      case 10 : return s * count > 500;
      case 20 : return s * count > 200;
      case 40 : return s * count > 100;
      case 100 : return s * count > 50;
      default : return s * count > 20;
      }
  }

private:
  Options* mOptions;
  bool mIsRead2;
  long mReads;
  int mEvaluatedSeqLen;
  /*
    why we use 8 here?
    map A/T/C/G/N to 0~7 by their ASCII % 8:
    'A' % 8 = 1
    'T' % 8 = 4
    'C' % 8 = 3
    'G' % 8 = 7
    'N' % 8 = 6
  */
  long[][8] mCycleQ30Bases; // might be long[][8]
  long[][8] mCycleQ20Bases;
  long[][8] mCycleBaseContents;
  long[][8] mCycleBaseQual;
  long[] mCycleTotalBase;
  long[] mCycleTotalQual;
  long[] mKmer;

  double[][string] mQualityCurves; // this might be an double[][string]
  double[][string] mContentCurves; // this might be an double[][string]
  long[string] mOverRepSeq;
  long[][string] mOverRepSeqDist; // this might be long[][string]

  size_t mCycles;
  size_t mBufLen;
  long mBases;
  long[8] mQ20Bases;
  long[8] mQ30Bases;
  long[8] mBaseContents;
  long mQ20Total;
  long mQ30Total;
  bool summarized;
  long mKmerMax;
  long mKmerMin;
  size_t mKmerBufLen;
}

unittest
{

}

unittest
{
    int r = 64, g = 121, b = 13;

    string kmer = "AGCTTGCGATC";
    long val = 643385;
    int prop = 50;

    char[] ss;
    ss.length = 255;

    int slen = sprintf(&ss[0], "<td style='background:#%.2X%.2X%.2X 'title='%s: %d\n %d times as mean value'> </td>",
                       r,  g,  b,           &kmer[0],val, prop);
    printf("%s -- ,len = %d ", &ss[0], slen);

}
