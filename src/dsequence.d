// Written in the D programming language.
import std.stdio;
///
class Sequence
{
  string mStr;

  this()
  { }

  this(string seq)
  {
    mStr = seq;
  }

  void print()
  {
    writef("%s", mStr);
  }

  size_t length()
  {
    return mStr.length;
  }
  /// Reverses complement
  Sequence reverseComplement(){
    char[] str;
    str.length = mStr.length;
    immutable len = mStr.length;
    foreach (size_t c ; 0 .. length)
      {
        char base = mStr[c]; 
        switch (base)
          {
          case 'A':
          case 'a':
            str[len-c-1] = 'T';
            break;
          case 'T':
          case 't':
            str[len-c-1] = 'A';
            break;
          case 'C':
          case 'c':
            str[len-c-1] = 'G';
            break;
          case 'G':
          case 'g':
            str[len-c-1] = 'C';
            break;
          default:
            str[len-c-1] = 'N';
          }
      }
    return new Sequence(cast(immutable) str);
  }

  /**
   * Operator overload
   *
   */
  Sequence opUnary(string s) () if (s == "~"){
    return reverseComplement();
  }
}

unittest
{
  auto s = new Sequence("AAAATTTTCCCCGGGG");
  Sequence rc = ~s;
  // auto rc = s.reverseComplement;
  assert(s.mStr == "AAAATTTTCCCCGGGG", "Failed in reverseComplement() expect AAAATTTTCCCCGGGG, but get " ~ s.mStr);
  assert(rc.mStr == "CCCCGGGGAAAATTTT", "Failed in reverseComplement() expect CCCCGGGGAAAATTTT, but get " ~ rc.mStr);
}
