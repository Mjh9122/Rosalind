class Fasta_dna:
    def __init__(self, lines):
        assert lines[0][0] == '>'
        self.id = lines[0][1:].strip()
        self.dna_string = ""
        for line in lines[1:]:
            self.dna_string += line.strip()

    def __repr__(self):
        return self.id + " " + self.dna_string
    
    def count_acgt(self):
        acgt = [0, 0, 0, 0]
        for char in self.dna_string:
            if char in {"A", 'a'}:
                acgt[0] += 1
            elif char in {'C', 'c'}:
                acgt[1] += 1
            elif char in {'G', 'g'}:
                acgt[2] += 1
            elif char in {'T', 't'}:
                acgt[3] += 1
        return acgt[0], acgt[1], acgt[2], acgt[3]
    
    def to_rna(self):
        return self.dna_string.replace('T', 'U')
    
    def reverse_compliment(self):
        rev = self.dna_string[::-1]
        molecule_map = {"A":"T", "T":"A","C":"G","G":"C"}
        rev_comp = ""
        for char in rev:
            rev_comp += molecule_map[char]
        return rev_comp

    def gc_content(self):
        tot = len(self.dna_string)
        cgs = 0
        for char in self.dna_string:
            if char in ("C", 'G'):
                cgs += 1
        return (cgs/tot) * 100
    
    def motif_locations(self, substring):
        locations = []
        length = len(substring)
        for index in range(len(self.dna_string)-length):
            if self.dna_string[index:index+length] == substring:
                locations.append(index + 1)
        return locations
    
    def hamming_dist(self, comparison_string: str) -> int:
        if len(self.dna_string) != len(comparison_string):
            raise Exception("DNA lengths must be the same for hamming distance")
        else:
            diff = 0
            for index, char in enumerate(self.dna_string):
                if char != comparison_string[index]:
                    diff += 1
        return diff


class Fasta_Batch_Load:
    def __init__(self, filepath):
        self.dna_list = []
        with open(filepath) as f:
            lines = f.readlines()
            start, end = 0, 0
            for index, line in enumerate(lines):
                if line[0] == '>':
                    end = index
                    if start != end:
                        self.dna_list.append(Fasta_dna(lines[start:end]))
                    start = index
            end = len(lines)
            self.dna_list.append(Fasta_dna(lines[start:end]))

    
    def __repr__(self) -> str:
        return str([dna.id for dna in self.dna_list])


