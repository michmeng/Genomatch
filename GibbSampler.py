import random

def hamming_dist(str1, str2):
    mismatch = 0
    #should prob check if len(str1) == len(str2)
    for i in range(len(str1)):
        c1 = str1[i]
        c2 = str2[i]
        if c1 != c2:
            mismatch += 1
    return mismatch

def most_freq_words(text, k):
    all_pat = {}
    for i in (range(len(text) - k + 1)):
        pattern = text[i:(i+k)]
        if pattern in all_pat:
            all_pat[pattern] = all_pat.get(pattern) + 1
        else:
            all_pat[pattern] = 1
    result = []
    max = 0
    for i in all_pat:
        if all_pat[i] > max:
            max = all_pat[i]
    for i in all_pat:
        if all_pat[i] == max:
            return i

def score(motifs):
    #find consensus string
    con = ""
    for i in range(len(motifs[0])):
        l = []
        for each in motifs:
            l.append(each[i])
        line = ""
        for each in l:
            line = line + each
        res = most_freq_words(line, 1)
        con = con + res
    #sum hamming dist between consensus string and each of motif
    sum = 0
    for each in motifs:
        diff = hamming_dist(each, con)
        sum += diff
    return sum


def profile_most_prob(text, k, profile):
    all_kmer = []
    for i in range(len(text) - k + 1):
        check = text[i:(i+k)]
        all_kmer.append(check)
    probs = {}
    for each in all_kmer:
        prob = 1.0
        index = 0
        for c in each:
            if c == "A":
                prob = float(prob * profile[index][0])
            if c == "C":
                prob = float(prob * profile[index][1])
            if c == "G":
                prob =  float(prob * profile[index][2])
            if c == "T":
                prob = float(prob * profile[index][3])
            index += 1
        probs[each] = float(prob)
    result = []
    max = 0
    for each in probs:
        if probs[each] > max:
            max = probs[each]
    for each in probs:
        if probs[each] == max:
            result.append(each)
    # if more than one result, choose the first one to appear
    min = len(probs) - 1
    if len(result) == 1:
        r = result[0]
    if len(result) > 1:
        for each in result:
            index = all_kmer.index(each)
            if index < min:
                min = index
        r = all_kmer[min]
    return r

def random_with_weight(lst):
    total = 0
    lst.insert(0,0)
    for each in lst:
        total += each
    #normalize if total does not equal one
    if total != 1:
        for i in range(len(lst)):
            if i == 0:
                lst[i] = float(lst[i]/total)
            else:
                 lst[i] = float((lst[i]/total)+lst[i-1])
    num = random.random()
    #find which index to return
    index = -1
    for i in range(len(lst)-1):
        if num >= lst[i] and num < lst[i+1]:
            index = i
    return index

def make_profile(motifs):
    columns = [''.join(seq) for seq in zip(*motifs)]
    return [[float(col.count(nuc) + 1) / float(len(col)) for nuc in 'ACGT'] for col in columns]

def generate_random_profile(profile, text, k):
    k_mers = []
    for i in range(0, len(text)-k+1):
        check = text[i:i+k]
        k_mers.append(check)
    probs = {}
    for each in k_mers:
        prob = 1.0
        index = 0
        for c in each:
            if c == "A":
                prob = float(prob * profile[index][0])
            if c == "C":
                prob = float(prob * profile[index][1])
            if c == "G":
                prob =  float(prob * profile[index][2])
            if c == "T":
                prob = float(prob * profile[index][3])
            index += 1
        probs[each] = float(prob)
    lst = []
    results = []
    for each in probs:
        results.append(each)
        lst.append(probs.get(each))
    index = random_with_weight(lst)
    return results[index]

def gibb_sampler(dna, k, t, N):
    # randomly select kmers to make motif, one from each string in DNA
    motifs = []
    for each in dna:
        patterns = []
        for i in range(len(each) - k + 1):
            patterns.append(each[i:(i+k)])
        index = random.randrange(0, len(patterns), 1)
        motifs.append(patterns[index])
    best_motif = motifs[:]

    for j in range(1, N+1):
        rand = random.randrange(0, t)
        curr = motifs.copy()
        del curr[rand]
        profile = make_profile(curr)
        motifs[rand] = generate_random_profile(profile, dna[rand], k)
        if score(motifs) < score(best_motif):
            best_motif = motifs
    return best_motif

def solver(Dna, k, t, N):
    i = 0

    last_motif = gibb_sampler(Dna, k, t, N)

    while(i<21):
        best_motif = gibb_sampler(Dna, k, t, N)
        if score(best_motif) < score(last_motif):
            last_motif = best_motif
        i+=1

    result = ""
    for each in last_motif:
        result = result + "\n" + each
    return result

    return last_motif

with open('/Users/michmeng/Documents/Genomatch/hirmed1.gDNA.longreads.fa') as f:
    lines = f.readlines()

k = 20
t = 15
N = 10

result = [i for i in lines if ">" not in i and len(i) > k]

# result = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
# "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
# "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
# "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
# "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

print (solver(result,k, t, N))