import std.stdio;
import etc.c.zlib;
import core.stdc.string: memset, strlen;
import core.stdc.stdio;
import util;
import std.conv : to;
//import utils.zio;
import read;
import common;

class FastqReader{
private:
  string mFilename;
  gzFile mZipFile;
  FILE* mFile;
  bool mZipped;
  bool mHasQuality;
  bool mPhred64;
  size_t _bytesRead;
  size_t _bytesTotal;
public:
  FastqReader* mLeft;
  FastqReader* mRight;
  this(string filename, bool hasQuality = true, bool phred64=false){
    mFilename = filename;
    mZipFile = null;
    mZipped = false;
    mPhred64 = phred64;
    mHasQuality = hasQuality;
    init();
  }
  ~this(){
    close();
  }
  bool isZipped(){
    return mZipped;
  }
  
  void getBytes(ref size_t bytesRead, ref size_t bytesTotal){
    if(mZipped) {
      bytesRead = gzoffset(mZipFile);
    } else {
      bytesRead = _bytesRead;
    }
    
    bytesTotal = _bytesTotal;
  }
  
  /// Note:
  ///this function is not thread-safe
  ///do not call read() of a same FastqReader object from different threads concurrently
  Read read(){
    if (mZipped){
      if (mZipFile == null || gzeof(mZipFile))
        return null;
    }
    assert(mZipFile || mFile !is null, "File: "~ mFilename ~ " could not be opend or something ... mZipped: ");

    if(!mFile || mFile.feof()) {
      return null;
    }

    string name = getLine();
    string sequence = getLine();
    string strand = getLine();

    if(name == "" || sequence == "" || strand == "")
      return null;

    // WAR for FQ with no quality
    if (!mHasQuality){
      import std.range: repeat;
		  import std.conv: to;
		  string quality = 'K'.repeat(sequence.length).to!string;
      return new Read(name, sequence, strand, quality, mPhred64);
    }
    else {
      string quality = getLine();
      if(quality == "")
        return null;
      return new Read(name, sequence, strand, quality, mPhred64);
    }

    return null;
  }
  bool eof(){
    if (mZipped) {
      return gzeof(mZipFile) != 0;
    } else {
      return mFile.feof() != 0;
    }
  }

  static bool isZipFastq(string filename){
    if (ends_with(filename, ".fastq.gz"))
      return true;
    else if (ends_with(filename, ".fq.gz"))
      return true;
    else if (ends_with(filename, ".fasta.gz"))
      return true;
    else if (ends_with(filename, ".fa.gz"))
      return true;
    else
      return false;
  }
  static bool isFastq(string filename){
    if (ends_with(filename, ".fastq"))
      return true;
    else if (ends_with(filename, ".fq"))
      return true;
    else if (ends_with(filename, ".fasta"))
      return true;
    else if (ends_with(filename, ".fa"))
      return true;
    else
      return false;
  }

private:
  void init() {
    if (ends_with(mFilename, ".gz")){
      mZipFile = gzopen(&mFilename[0], "rb".ptr);
      mZipped = true;
      gzrewind(mZipFile);
    }
    if (isFastq(mFilename)){
      mFile = fopen(&mFilename[0], "rb".ptr);
      assert(mFile, "fopen failed");
      
      fseek(mFile, 0, SEEK_END);
      _bytesTotal = ftell(mFile);
      fseek(mFile,0, SEEK_SET);
      
      mZipped = false;
      
    }
  }
  void close(){
    if (mZipped){
      if (mZipFile){
        gzclose(mZipFile);
        mZipFile = null;
      }
    }
    else {
      if (mFile){
        mFile.fclose();
        mFile = null;
      }
    }
  }
  string getLine(){
      static immutable size_t maxLineLength = 2048;
      char[maxLineLength] line = 0;
      char* linep = &line[0];
      printf("Start _bytesRead: %d, strlen(linep) %d \"%s\"\n",_bytesRead, strlen(linep), linep);
      if(mZipped) {
        char* buf = null;
          memset(linep, 0, maxLineLength);
          buf = cast(char*) gzgets(mZipFile, linep, maxLineLength);
    
          // EOF or error, return an empty string
          if(!buf)
            return "";
    
          // optimize for short reads
          // reached end of line
          if(line[maxLineLength-2] == '\0' || line[maxLineLength-2] == '\n') {
            clearLineBreaks(linep);
  
            return line[0 .. strlen(linep)].to!string;
          } else {
            string s;
            while(buf) {
              memset(linep, 0, maxLineLength);
              buf = cast (char*) gzgets(mZipFile, linep, maxLineLength);
              //eof or error
              if(!buf)
                break;
              //reached end of line
              if(line[maxLineLength-2] == '\0' || line[maxLineLength-2] == '\n') {
              clearLineBreaks(buf);
              s ~= buf[0 .. strlen(buf)];
              break;
            } else {
              s ~= buf[0 .. strlen(buf)];
            }
          }
          return s;
        }
      }
      else {
        char* lptr = linep;
        // Because we are not using stream, we need to keep track of number of bytes that have been read.
        printf("_bytesRead: %d, strlen(linep) %d \"%s\"\n",_bytesRead, strlen(linep), linep);
        _bytesRead += getline(&lptr, cast(size_t*)&maxLineLength, mFile);
        // optimize for short reads
        // reached end of line
        if(mFile && mFile.feof()) {
          clearLineBreaks(linep);
          printf("_bytesRead: %d, strlen(linep) %d \"%s\"\n",_bytesRead, strlen(linep), linep);
          return line[0 .. strlen(linep)].to!string;
        } else {
          string s = line[0 .. strlen(linep)].to!string;
          while(true) {
            if(!mFile || mFile.feof()) 
              break;
            memset(linep, 0, maxLineLength);
            getline(&lptr, cast(size_t*)&maxLineLength, mFile);
            
            if(mFile && mFile.feof()) {
              clearLineBreaks(linep);
              s ~= line[0 .. strlen(linep)];
              break;
            } else {
              // in case of some error happened, break it
              if(line[0] == '\0'){
                break;
              }
              s ~= line[0 .. strlen(linep)];
            }
          }
          return s;
        }
      }
  
      return "";
    }
  size_t clearLineBreaks (char* line) {
    // trim \n, \r or \r\n in the tail
    ulong readed = strlen(line);
    if(readed >= 2) {
      if(line[readed-2] == '\n' || line[readed-1] == '\r'){
        line[readed-2] = '\0';
        readed -= 1;
      }
        else  if(line[readed-2] == '\r'|| line[readed-2] == '\n') {
          line[readed-2] = '\0';
          readed -= 2;
      }
    }
    return readed;
  }
}

class FastqReaderPair{
  FastqReader mLeft;
  FastqReader mRight;

  this(FastqReader left, FastqReader right){
    mLeft = left;
    mRight = right;
  }
  this(string leftName, string rightName, bool hasQuality = true, bool phred64 = false){
    mLeft = new FastqReader(leftName, hasQuality, phred64);
    mRight = new FastqReader(rightName, hasQuality, phred64);
  }
  ~this(){
    if(mLeft){
      mLeft.destroy;
      mLeft = null;
    }
    if(mRight){
      mRight.destroy;
      mRight = null;
    }
  }
  ReadPair read(){
    Read l = mLeft.read();
    Read r = mRight.read();
    if(!l || !r){
      return null;
    } else {
      return new ReadPair(l, r);
    }
  }

}

unittest{
  bool test(){
   auto reader1 = new FastqReader("testdata/R1.fq");
   auto reader2 = new FastqReader("testdata/R1.fq.gz");
    Read r1 = null;
    Read r2 = null;
    while(true){
      r1=reader1.read();
      r2=reader2.read();
      if(r1 is null || r2 is null)
        break;
      //printf("r1.mSeq.mStr%*.s:\nr2.mSeq.mStr%*.s", r1.mSeq.mStr.length, r1.mSeq.mStr.ptr, r2.mSeq.mStr.length, r2.mSeq.mStr.ptr);
      if (r1.mSeq.mStr != r2.mSeq.mStr){
        return false;
      }
      r1.destroy;
      r2.destroy;
	}
	return true;
  }
  assert(test());

}
