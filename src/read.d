import dsequence;
import util;
import std.stdio: File, writeln;
import std.algorithm.comparison: min, max;
import std.algorithm.mutation: reverse;
import std.string;

class Read{
  string mName;
  Sequence mSeq;
  string mStrand;
  string mQuality;
  bool mHasQuality;

  this(string name, string seq, string strand, string quality, bool phred64=false){
    mName = name;
    mSeq = new Sequence(seq);
    mStrand = strand;
    mQuality = quality;
    mHasQuality = true;
    if(phred64)
      convertPhred64To33();

  }
  this(string name, Sequence seq, string strand, string quality, bool phred64=false){
    mName = name;
    mSeq = seq;
    mStrand = strand;
    mQuality = quality;
    mHasQuality = true;
    if(phred64)
      convertPhred64To33();

  }
  this(string name, string seq, string strand){
    mName = name;
    mSeq = new Sequence(seq);
    mStrand = strand;
    mHasQuality = false;

  }
  this(string name, Sequence seq, string strand){
    mName = name;
    mSeq = seq;
    mStrand = strand;
    mHasQuality = false;

  }
  this(ref Read r){

    mName = r.mName;
    mSeq = r.mSeq;
    mStrand = r.mStrand;
    mQuality = r.mQuality;
    mHasQuality = r.mHasQuality;
  }

  void print(){
    writeln( mName );
    writeln( mSeq.mStr );
    writeln( mStrand );
    if(mHasQuality)
      writeln( mQuality );

  }
  void printFile(File file){
    file.writeln( mName );
    file.writeln( mSeq.mStr);
    file.writeln( mStrand );
    if(mHasQuality)
      file.writeln( mQuality );

  }
  Read reverseComplement(){
    Sequence seq = ~mSeq;
    dchar[] qual = cast(dchar[])mQuality.dup;
    reverse(qual);
    string strand = (mStrand=="+") ? "-" : "+";
    Read newRead = new Read(mName, seq, strand, cast(string)qual);
    return newRead;

  }
  string firstIndex(){
    size_t len = this.mName.length;
    size_t end = len;
    if(len<5)
      return "";
    for(size_t i=len-3;i>=0;i--){
      if(mName[i]=='+')
        end = i-1;
      if(mName[i]==':'){
        return mName[i+1 .. end];
      }
    }
    return "";

  }
  string lastIndex(){
    size_t len = mName.length;
    if(len<5)
      return "";
    for(size_t i=len-3;i>=0;i--){
      if(mName[i]==':' || mName[i]=='+'){
        return mName[i+1 .. len];
      }
    }
    return "";

  }
  // default is Q20
  size_t lowQualCount(int qual=20){
    size_t count = 0;
    for(int q=0;q<mQuality.length;q++){
      if(mQuality[q] < qual + 33)
        count++;
    }
    return count;

  }
  size_t length(){
    return mSeq.length;
  }
  override string toString(){
    return mName ~ "\n" ~ mSeq.mStr ~ "\n" ~ mStrand ~ "\n" ~ mQuality ~ "\n";
  }
  void resize(int len){
    if(len > this.length() || len<0)
      return ;
    mSeq.mStr = mSeq.mStr[0..len];
    mQuality = mQuality[0..len];
  }
  void trimFront(int len){
    len = min(this.length()-1, len);
    mSeq.mStr = mSeq.mStr[len .. mSeq.mStr.length - len];
    mQuality = mQuality[len .. mQuality.length - len];

  }

  void convertPhred64To33(){
    dchar[] ret;
    for(size_t i=0; i<mQuality.length; i++) {
      ret[i] = max(33, mQuality[i] - (64-33));
    }
    this.mQuality = cast(string)ret;
  }
}

unittest{
  bool test(){
    Read r;
	r = new Read("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
           "CTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTCCTTAGGAGGACATTTTTTACATGAAATTATTAACCTAAATAGAGTTGATC",
           "+",
           "AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEEEAEEEEEEEEEEEEEEEAEEE/EEEEEEEEEEAAEAEAAEEEAEEAA");
	string idx = r.lastIndex();
	return idx == "GGTCCCGA";

  }

  assert(true == test());

}

class ReadPair{
  Read mLeft;
  Read mRight;

  this(ref Read left, ref Read right){
    mLeft = left;
	mRight = right;
  }
  ~this(){
    if(mLeft){
      delete mLeft;
      mLeft = null;
	}
	if(mRight){
      delete mRight;
      mRight = null;
	}

  }

    // merge a pair, without consideration of seq error caused false INDEL
    Read fastMerge(){
	Read rcRight = mRight.reverseComplement();
	size_t len1 = mLeft.length();
	size_t len2 = rcRight.length();
	// use the pointer directly for speed
	const char* str1 = mLeft.mSeq.mStr.ptr;
	const char* str2 = rcRight.mSeq.mStr.ptr;
	const char* qual1 = mLeft.mQuality.ptr;
	const char* qual2 = rcRight.mQuality.ptr;

	// we require at least 30 bp overlapping to merge a pair
	const size_t MIN_OVERLAP = 30;
	bool overlapped = false;
	size_t olen = MIN_OVERLAP;
	size_t diff = 0;
	// the diff count for 1 high qual + 1 low qual
	size_t lowQualDiff = 0;

	while(olen <= min(len1, len2)){
		diff = 0;
		lowQualDiff = 0;
		bool ok = true;
		size_t offset = len1 - olen;
		for(size_t i=0;i<olen;i++){
			if(str1[offset+i] != str2[i]){
				diff++;
				// one is >= Q30 and the other is <= Q15
				if((qual1[offset+i]>='?' && qual2[i]<='0') || (qual1[offset+i]<='0' && qual2[i]>='?')){
					lowQualDiff++;
				}
				// we disallow high quality diff, and only allow up to 3 low qual diff
				if(diff>lowQualDiff || lowQualDiff>=3){
					ok = false;
					break;
				}
			}
		}
		if(ok){
			overlapped = true;
			break;
		}
		olen++;
	}

	if(overlapped){
		size_t offset = len1 - olen;
		string mergedName =  mLeft.mName ~ " merged offset:" ~
          format("%s, overlap: %s, diff: %s", offset, olen, diff);
		auto mergedSeq = cast(dchar[])(mLeft.mSeq.mStr[0..offset] ~ rcRight.mSeq.mStr);
		auto mergedQual = cast(dchar[])(mLeft.mQuality[0..offset] ~ rcRight.mQuality);
		// quality adjuction and correction for low qual diff
		for(size_t i=0;i<olen;i++){
			if(str1[offset+i] != str2[i]){
				if(qual1[offset+i]>='?' && qual2[i]<='0'){
					mergedSeq[offset+i] = str1[offset+i];
					mergedQual[offset+i] = qual1[offset+i];
				} else {
					mergedSeq[offset+i] = str2[i];
					mergedQual[offset+i] = qual2[i];
				}
			} else {
				// add the quality of the pair to make a high qual
              mergedQual[offset+i] =  cast(dchar)(qual1[offset+i] + qual2[i] - 33);
			}
		}
		delete rcRight;
		return new Read(cast(string)mergedName, cast(string)mergedSeq, "+", cast(string)mergedQual);
	}

	delete rcRight;
	return null;
}
}

unittest{
	Read left = new Read("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
                          "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAG",
                          "+",
                          "AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
	Read right = new Read("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
                           "AAAAAACTACACCATAGAATGACTATGAGTCTCATAAGAATGCACTCAACTAGTCATCACTCCTGTGTTTTCATAAGAAAAAACAGTGTTAGAGTCCAAGAG",
                           "+",
                           "AAAAA6EEEEE/EEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");

	auto pair = new ReadPair(left, right);
	Read merged = pair.fastMerge();
	// assert(null !is merged);
	// assert(merged.mSeq.mStr == "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTTT");
}
