import std.string;
immutable string[string] KnownAdapters;
static this()
{
  KnownAdapters = [
                   "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA": "Illumina TruSeq Adapter Read 1",
                   "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT": "Illumina TruSeq Adapter Read 2"
                   ];


}

unittest{
  assert("Illumina TruSeq Adapter Read 1" == KnownAdapters["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"]);
  assert("Illumina TruSeq Adapter Read 2" == KnownAdapters["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"]);
}
