// Written in the D programming language.
import std.stdio;
import dunittest;

extern (C++) int cppmain();

enum FASTP_VER = "0.13.2";

int main(string[] args)
{
    // display version info if no argument is given
    if (args.length == 1) {
        writefln("fastp: an ultra-fast all-in-one FASTQ preprocessor\nversion %s", FASTP_VER);
    }
    if (args.length == 2 && args[1] == "test"){
        UnitTest tester;
        tester.run();
        return 0;
    }
    if (args.length == 2 && (args[1] == "-v" || args[1] == "--version")){
        writefln("fastp: an ultra-fast all-in-one FASTQ preprocessor\nversion %s", FASTP_VER);
        return 0;
    }

    return cppmain();
}
