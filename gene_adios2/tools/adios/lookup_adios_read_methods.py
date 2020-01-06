#!/usr/bin/env python


import re
import argparse
import os


def FindInts(afile):
    f = open(afile, "r")
    text = f.read()
    f.close()

    test = re.compile("ADIOS_READ_METHOD_.*=\s*(\d+)")
    matches = test.findall(text)

    ints = []
    for match in matches:
        ints.append(int(match))

    i_min = min(ints)
    i_max = max(ints)
    return i_min, i_max


def Replace(mstr, text, index, join="\n"):
    #test = re.compile("(#define)\s*({0}_READ_METHOD)(.*)".format(mstr))
    #rep = "#define {0}_READ_METHOD {1}".format(mstr, index)
    #outtext = test.sub(rep, text)

    rep = "#define {0}_READ_METHOD {1}".format(mstr, index)
    outtext = "{0}{2}{1}".format(text, rep, join)
    return outtext

    return outtext


if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--adios",   help="ADIOS code",   required=True)
    parser.add_argument("-o", "--outfile", help="output file",  required=True)
    args = parser.parse_args()

    rfile = os.path.join(args.adios, "include", "adios_read_v2.h")
    min_index, max_index = FindInts(rfile)

    outtext = ""
    outtext = Replace("MIN", outtext, min_index, join="")
    outtext = Replace("MAX", outtext, max_index, join="\n")

    f = open(args.outfile, "w")
    f.write(outtext)
    f.close()
