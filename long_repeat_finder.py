import string
import sys
import cs190FileUtil # Utility to read Fasta files

def findLongestRepeat(text):
    max = 1 # Our current goal
    maxPos = -1
    maxDup = -1

     # Start at the begining, and continue for each spot
    for pos in range(len(text)):
        # Look for a match to the string we are sitting on
        dup = text.find(text[pos:pos+max], pos+1, len(text))

        # We have a match: can we extend it?
        while (dup > 0):
            maxPos = pos
            maxDup = dup
            max = max + 1 # Now look for a longer match

            # Can we find a longer match?
            dup = text.find(text[pos:pos+max], dup, len(text))

    return [maxPos, maxDup, max-1]

if (len(sys.argv) != 2):
    print ("Usage: python", sys.argv[0], "<filename>")
else:
    text = cs190FileUtil.readFastaFile(sys.argv[1])
    [pos, dup, ln] = findLongestRepeat(text)
    print ("Found duplicate of length", ln)
    print (pos, text[pos:pos+ln])
    print (dup, text[dup:dup+ln])