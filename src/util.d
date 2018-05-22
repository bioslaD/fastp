import std.stdio;
import std.string;
import std.array;
import std.algorithm;
import std.algorithm.iteration : map;
import std.algorithm.mutation : copy;
import std.array : appender;
import core.stdc.stdlib: exit;
import std.ascii: toLower, isAlpha;
import std.path: dirName, isDir;
import std.file: exists;

char complement(char base) {
    switch(base){
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        default:
            return 'N';
    }
}

bool starts_with( const ref string value, const ref string starting) {
    if (starting.length > value.length) return false;
    return startsWith(value, starting) == 1;
}

bool ends_with( const ref string value, const ref string ending)
{
	if (ending.length > value.length) return false;

	return  endsWith(value, ending) == 1;
}

string trim(const ref string str)
{
  return str.strip;
}

/**
 * This split function is probably never used.
 */

int split(const ref string str, ref string[] ret_, string sep = ",")
{
    if (str.empty())
    {
        return 0;
    }

    string tmp;
    size_t pos_begin = str.indexOfNeither(sep);
    size_t comma_pos = 0;

    while (pos_begin != -1)
    {
      comma_pos = str.indexOf(sep, pos_begin);
        if (comma_pos != -1)
        {
          tmp = str[pos_begin .. comma_pos];
            pos_begin = comma_pos + sep.length;
        }
        else
        {
            tmp = str[pos_begin..$];
            pos_begin = comma_pos;
        }

        ret_ ~= tmp;
        delete(tmp);
    }
    return 0;
}

/**
 * replace all occurrences of `src` with `dest` that exists in `str`.
 *
 * Returns: a new string
 */

string replace(const ref string str, const ref string src, const ref string dest)
{
    string ret;

    size_t pos_begin = 0;
    size_t pos       = str.indexOf(src);
    while (pos != -1)
    {
        ret ~= str[pos_begin .. pos];
        ret ~= dest;
        pos_begin = pos + src.length;
        pos       = str.indexOf(src, pos_begin);
    }
    if (pos_begin < str.length)
    {
        ret ~= str[pos_begin..$];
    }
    return ret;
}

/**
 * This function uses `char` (8bit) instead of `dchar` (32bit)
 * So it expects input file to be ascii.
 *
 */

string reverse(const ref string str) {
    char[] ret = new char[str.length];
    for(int pos=0; pos<str.length; pos++) {
        ret[pos] = str[str.length - pos - 1];
    }
    return cast(string)ret;
}

void check_file_valid(const ref string s) {
    if(!s.exists){
      stderr.writeln("ERROR: file '", s, "' doesn't exist, quit now");
        exit(-1);
    }
    if(s.isDir){
      stderr.writeln("ERROR: '" , s, "' is a folder, not a file, quit now");
        exit(-1);
    }
}
//
void check_file_writable(const ref string s) {
    string dir = dirName(s);
    if(!dir.exists) {
       stderr.writeln( "ERROR: '", dir, " doesn't exist. Create this folder and run this command again.");
        exit(-1);
    }
    if(s.isDir){
        stderr.writeln("ERROR: '", s, "' is not a writable file, quit now");
        exit(-1);
    }
}
//
/// Remove non alphabetic characters from a string
string str_keep_alpha(const  ref string s)
{
     string new_str;
    for( size_t it =0; it < s.length; it++) {
        if(  isAlpha(s[it]) ) {
            new_str ~= s[it];
        }
    }
    return new_str;
}

/// Remove invalid sequence characters from a string
string str_keep_valid_sequence(const ref string s)
{
     string new_str;
    for( size_t it =0; it < s.length; it++) {
        if(  isAlpha(s[it]) || s[it] == '-' || s[it] == '*' ) {
            new_str ~= s[it];
        }
    }
    return new_str;
}

int find_with_right_pos(const ref string str, const ref string pattern, int start=0) {
    int pos = cast(int)str.indexOf(pattern, start);
    if (pos < 0)
        return -1;
    else
        return pos + cast(int)pattern.length;
}

void changeCase(alias F)(ref string s){
  // auto abuf = appender!(char[])();
  // s.map!F.copy(abuf);
  // s = cast(string)abuf.data;
  string ret;
  foreach(ref c; s){ ret ~= F(c);}
  s = ret;

}

void str2upper(ref string s){
     changeCase!toUpper(s);
}

/// This function expects ASCII input
void str2lower(ref string s){
    changeCase!toLower(s);
 
}

char num2qual(int num) {
    if(num > 126 - 33)
        num = 126 - 33;
    if(num < 0)
        num = 0;

    char c = cast(char)(num + 33);
    return c;
}

void error_exit(const string msg) {
  writeln("ERROR: ", msg );
    exit(-1);
}

unittest{
  writeln("Test complement");
  string cpls = "AtGC";
  // assert("TACG" == cpls.map!(c => complement(cast(char)c)).data);
  assert('T' == complement('a'));

  string a = "abc";
  string b = "ab";
  string c = "b";
  writeln("Test starts_with");
  assert(starts_with(a, b));
  assert(!starts_with(a, c));

  writeln("Test ends_with");
  string d = "c";
  string e = "d";
  assert(ends_with(a, d));
  assert(!ends_with(a, e));

  writeln("Test split");
  string spl = "|abc|def|gh";
  string[] ret;
  split(spl, ret, "|");
  assert(equal(ret, ["abc", "def", "gh"]));
  writeln(ret);
  delete(ret);

  writeln("Test replace");
  auto srcstr = "Test replacement test";
  auto src = "est";
  auto dest = "esting";
  auto rpr = replace(srcstr, src, dest);
  assert(rpr == "Testing replacement testing");

  writeln("Test no replacement");
  assert(replace(src, dest, dest) == "est");

  assert(num2qual(83) == 't', "Expected output is 't', but got " ~ num2qual(83));
  assert(num2qual(40) == 'I', "Expected output is 'I', but got " ~ num2qual(40));
  assert(num2qual(1000) == '~', "Expected output is '~', but got " ~ num2qual(1000));

  string tlower = "Test LOWER";
  string tlower2 = "Test LOWER2";
  writeln(tlower, " Before str2lower");
  str2lower(tlower);
  changeCase!toLower(tlower2);
  writeln(tlower.map!(c =>std.ascii.toLower(c)), " After str2lower");
  writeln(tlower2, " After changeCase");

  writeln("Test reverse");
  string s = "";
  assert("" == reverse(s));
  s = "abc";
  assert("cba" == reverse(s));
  s = "aBc";
  assert("cBa" == reverse(s));

  // This fails because the function does not know how to handle unicode yet.
  s = "ư";
  assert("ư" == reverse(s));

}

/// Benchmarking
unittest{
  import std.datetime.stopwatch: benchmark, Duration;
  import std.meta;
  import std.range: iota;
  import std.conv: to;
  string randomDna(int length) {
    import std.random;
    auto rnd = Random(unpredictableSeed);
    enum chars = ['a','C','g','T'];
    return iota(length).map!(a=>chars[uniform(0,4, rnd)]).array;
  }

  auto input = randomDna(3000);

  void bm_str2lower(){
    string tlower = input;
    str2lower(tlower);
  }
  alias m2lower = changeCase!toLower;
  void bm_changeCase(){
    string cc = input;
    m2lower(cc);
  }

  string prev = null;
  foreach(fn; AliasSeq!(bm_str2lower, bm_changeCase)){
    auto timing = benchmark!({fn();})(10_000);
    writeln((&fn).stringof[4..$], ": ", timing[0].to!Duration);
  }
}
