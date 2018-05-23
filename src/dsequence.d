// Written in the D programming language.
import std.stdio;

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

    Sequence reverseComplement(){
        char[] str = new char[mStr.length];
        for (size_t c = 0; c < mStr.length; c++)
        {
            char base = mStr[c];
            switch (base)
            {
                case 'A':
                case 'a':
                    str[mStr.length-c-1] = 'T';
                    break;
                case 'T':
                case 't':
                    str[mStr.length-c-1] = 'A';
                    break;
                case 'C':
                case 'c':
                    str[mStr.length-c-1] = 'G';
                    break;
                case 'G':
                case 'g':
                    str[mStr.length-c-1] = 'C';
                    break;
                default:
                    str[mStr.length-c-1] = 'N';
            }
        }
        return new Sequence(str.idup);
    }

/*
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
