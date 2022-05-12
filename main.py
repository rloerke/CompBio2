"""
Comp Bio Assignment 2
Written by Ray Loerke
Primitive Gene Finder and Predictor
"""


def main():
    predict("mixed_overlapping_1.fasta", True)


class Gene:
    """
    This class will be used to create an object for each gene found in the genome
    dna, rna, and protein will hold the corresponding sequences
    direction is used to indicate if a gene was found on the forward ("+") or reverse ("-") strand
    The stop and start attributes are used to keep track of the positioning of the dna
    """

    def __init__(self, dna, dna_start, dna_stop, direction):
        self.dna = dna
        self.rna = ""
        self.protein = []
        self.dna_start = dna_start
        self.dna_stop = dna_stop
        self.direction = direction


class Protein:
    """
    This class will be used to keep track of each protein translated from RNA
    sequence holds the list of proteins
    start and stop keep track of the positioning of the protein
    """

    def __init__(self, sequence):
        self.sequence = sequence
        self.start = 0
        self.stop = 0


def predict(filename, operons):
    """
    This function will scan a genome and find the genes it contains.
    It will then find the rna and protein made from these genes

    :param filename: The name of the fasta file containing the genome to be scanned
    :param operons: This bool determines if the genome should be searched for operons
    :return: Files will be created for the genes dna, rna, and protein
    """
    # A helper function is called to interpret the fasta file
    f_strand, g_name = read_fasta(filename)

    # A helper function is used to infer the reverse strand
    r_strand = reverse_strand(f_strand)
    r_strand = r_strand[::-1]

    # A helper function is called to search the forward and reverse strands
    gene_list = predict_helper(f_strand, True, operons)
    gene_list_r = predict_helper(r_strand, False, operons)

    # A helper function is called to write the information about the genes we have found to several files
    write_fafsa(gene_list + gene_list_r, g_name, filename)


def read_fasta(filename):
    """
    This functions reads in and formats a fasta file

    :param filename: The name of the fasta file
    :return: It outputs a string containing the sequence found in the file
    """
    # The fasta file is opened, read in as a string, then the file is closed
    fasta = open(filename, "r")
    fasta_text = fasta.read()
    fasta.close()

    # The string is split into a list of strings based on newline characters
    fasta_list = fasta_text.split('\n')

    # The name of the file is extracted from the description line
    g_name = fasta_list[0][2:]
    counter = len(g_name)

    # The length of the strand is removed from the file name
    while g_name[counter - 1].isdigit():
        g_name = g_name[:counter - 1]
        counter -= 1

    # This variable will hold the sequence
    sequence = ""

    # Loop through each line of the file
    for line in fasta_list:

        # If the line contains nucleotide bases add them to our working sequence
        if not line.startswith('>') and line != '':
            sequence += line

    # Return our the sequence
    return sequence, g_name


def reverse_strand(f_strand):
    """
    This function takes a DNA strand and finds its compliment

    :param f_strand: The DNA strand to be reversed
    :return: The complementary DNA strand
    """
    r_strand = ""

    # Each letter is replaced with its pair (A<->T) (C<->G)
    for letter in f_strand:
        if letter == 'A':
            r_strand += 'T'
        elif letter == 'T':
            r_strand += 'A'
        elif letter == 'C':
            r_strand += 'G'
        elif letter == 'G':
            r_strand += 'C'

    return r_strand


def gene_find(strand, gene_list, direc, length):
    """
    This function searches the DNA strand for a promoter and terminator.
    If both are found it creates a Gene object and recursively searches the rest of the strand for more genes

    :param strand: The DNA strand to be searched
    :param gene_list: The list of genes that have been found
    :param direc: Indicates if the strand is the forward (True) or reverse (False) strand
    :param length: The original length of the DNA strand (before it is cut upon recursive calls)
    :return: The list of Gene objects that have been found
    """
    # Variables to track the index where terminator start, if a hairpin loop exists, and the hairpin loop
    t_index = -1
    hairpin = False
    h_string = ""

    # The strand is searched for the Pribnow Box
    p_index = strand.find("TATAAT")

    # If the Pribnow Box was found move the index up to the start of the gene and search for the terminators poly-U tail
    if p_index != -1:
        p_index += 7
        t_index = strand.find("TTTTTTTTTT", p_index + 6)

    # If the terminator was found search before it for a GC-rich hairpin loop
    if t_index != -1:
        for it in range(t_index - 20, t_index - 12):
            if strand[it] == 'C' or strand[it] == 'G':
                h_string += strand[it]

        # If the first side of the hairpin loop is long enough check if it is complimented, if so we have found a gene
        if len(h_string) == 8:
            if strand[t_index - 8: t_index] == reverse_strand(h_string[::-1]):
                hairpin = True

    # The direction attribute of the Gene is set up
    if direc:
        d = "+"
    else:
        d = "-"

    num_genes = len(gene_list)

    # If this is the first gene we have found its start and stop are found by moving our indices forward
    if not gene_list:
        start_pos = p_index + 6
        stop_pos = t_index + 9

    # If we have found other genes the indices need to be moved forward based on where the last gene ended
    else:
        if direc:
            start_pos = p_index + 6 + gene_list[num_genes - 1].dna_stop
            stop_pos = t_index + 9 + gene_list[num_genes - 1].dna_stop

        # If we are on the reverse strand we need to reverse the dna_stop of the previous gene before adding it
        else:
            start_pos = p_index + 6 + length - gene_list[num_genes - 1].dna_stop
            stop_pos = t_index + 9 + length - gene_list[num_genes - 1].dna_stop

    # If we have found a gene create a Gene object
    # The strand direction is checked so that the correct positioning can be assigned
    if hairpin:
        if direc:
            g = Gene(strand[p_index + 6: t_index + 10], start_pos + 1, stop_pos + 1, d)
        else:
            g = Gene(strand[p_index + 6: t_index + 10], length - start_pos, length - stop_pos, d)

        # If there are more nucleotides in the strand make a recursive call to search for more genes
        if strand[t_index + 10:] != "":
            return gene_find(strand[t_index + 10:], gene_list + [g], direc, length)

        # If this is the last gene add it to our list and return that list
        else:
            gene_list += [g]

    # A new Gene list is created to hold only the unique genes
    new_gene_list = []
    for gene in gene_list:

        # If the gene is on the forward strand
        if direc:

            # If the gene wraps around or is from the first half of the doubled strand
            if gene.dna_start < length + 13 or gene.dna_stop < length:

                # If the terminator wraps around adjust the stop position
                if gene.dna_stop > length:
                    gene.dna_stop = gene.dna_stop - length

                # If the promoter wraps around adjust the start position
                if gene.dna_start > length:
                    gene.dna_start = gene.dna_start - length

                new_gene_list += [gene]

        # If the gene is on the reverse strand
        else:

            # If the gene wraps around or is from the first half of the double strand
            if gene.dna_start > -9 or gene.dna_stop > 0:

                # If the promoter wraps around adjust the stop position
                if gene.dna_stop < 0:
                    gene.dna_stop = gene.dna_stop + length

                # If the terminator wraps around adjust the start position
                if gene.dna_start < 0:
                    gene.dna_start = gene.dna_start + length

                new_gene_list += [gene]

    return new_gene_list


def translate(gene, operons):
    """
    This function translates a sequence of RNA into a list of proteins

    :param gene: The Gene containing the RNA to be translated
    :param operons: Determines if the gene should be searched for operons
    :return: Adds the protein list and protein positions to the Gene object, returns nothing
    """
    codon_list = []
    stop_index = -1
    rna = gene.rna

    # Check if this is the first protein
    if not gene.protein:

        # The RNA sequence is searched for the start codon
        start_index = rna.find("AUG")
    else:

        # Search for a start codon after the last proteins stop codon
        start_index = rna.find("AUG", gene.protein[len(gene.protein) - 1].stop)

    # If a start codon was found, search for one of the stop codons that is in-frame
    if start_index != -1:
        for it in range(start_index + 3, len(rna) - 3, 3):
            if rna[it:it + 3] == "UAA" or rna[it:it + 3] == "UAG" or rna[it:it + 3] == "UGA":
                stop_index = it
                break

    # If a stop codon was found, go through the RNA sequence and built a list of codons
    if stop_index != -1:
        for it in range(start_index, stop_index, 3):
            codon_list.append(rna[it:it + 3])

    protein_list = ""

    # This dictionary matches codons to their single letter protein classifications
    codon_dict = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
                  "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
                  "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
                  "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
                  "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                  "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                  "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                  "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                  "UAU": "Y", "UAC": "Y",
                  "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                  "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                  "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                  "UGU": "C", "UGC": "C", "UGG": "W",
                  "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                  "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                  "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

    # The codon list is put through the dictionary to get a protein list
    for codon in codon_list:
        protein_list += codon_dict[codon]

    # If a protein was found create a protein abject and add it to the genes protein list
    if start_index != -1 and stop_index != -1:
        p = Protein(protein_list)
        gene.protein += [p]

        # Set the start and end positions of the protein
        gene.protein[len(gene.protein) - 1].start = start_index + 1
        gene.protein[len(gene.protein) - 1].stop = stop_index

        # Perform a recursive call to look for operons if the option is enabled
        if operons:
            translate(gene, operons)


def write_fafsa(gene_list, gene_name, filename):
    """
    This function builds output files using the information in the list of Gene objects

    :param gene_list: List of Gene objects
    :param gene_name: The name of the gene pulled from the fasta file description line
    :param filename: The name of the fasta file containing the genome
    :return: Returns nothing, creates three output files:
    <filename>.g.fasta containing the DNA sequences of each gene
    <filename>.r.fasta containing the RNA sequences of each gene
    <filename>.p.fasta containing the Protein sequences of each gene
    """
    # The strings for each output file are created
    dna_output = ""
    rna_output = ""
    protein_output = ""

    # This protein counter will keep track of how many proteins we need to name
    p_counter = 0

    # Loop through each gene and add its info to each file
    for it in range(len(gene_list)):

        # A helper function is called to add the Genes information to each file
        # The corresponding suffix is passed in with an empty string indicating the sequence is DNA
        # The protein counter is also updated
        dna_add, p_counter = format_output(gene_list, gene_name, "", it, p_counter)
        dna_output += dna_add
        rna_add, p_counter = format_output(gene_list, gene_name, "r", it, p_counter)
        rna_output += rna_add
        protein_add, p_counter = format_output(gene_list, gene_name, "p", it, p_counter)
        protein_output += protein_add

    # The output files are created, the output strings are written into them, and the files are closed
    dna_file = open(filename[:len(filename) - 5] + "g.fasta", "w")
    rna_file = open(filename[:len(filename) - 5] + "r.fasta", "w")
    protein_file = open(filename[:len(filename) - 5] + "p.fasta", "w")
    dna_file.write(dna_output)
    rna_file.write(rna_output)
    protein_file.write(protein_output)
    dna_file.close()
    rna_file.close()
    protein_file.close()
    
    # Notify the used that gene finding is complete
    print("Gene Finding Completed!")


def format_output(gene_list, gene_name, suffix, it, p_counter):
    """
    This function creates the string to be written into the output files

    :param gene_list: The list of Genes
    :param gene_name: The name of the genome
    :param suffix: The suffix indicating what type of sequence we are working with
    :param it: An iterator indicating which gene in the list we are on
    :param p_counter: The number of proteins that need names
    :return: The string to be written into the output file
    """
    # Variables are created for each part of the description line and the sequence
    output = ""
    contents = ""
    direc = ""
    start = ""
    stop = ""

    # If we are working on the DNA output, all info is pulled from the Gene object
    if suffix == "":
        contents = gene_list[it].dna
        direc = gene_list[it].direction
        start = str(gene_list[it].dna_start)
        stop = str(gene_list[it].dna_stop)

    # If we are working on the RNA output
    if suffix == "r":
        contents = gene_list[it].rna

        # With RNA the strand direction is always positive and the position is based on the length of the sequence
        direc = "+"
        start = "1"
        stop = str(len(gene_list[it].rna))

    # If we are working on the protein output
    if suffix == "p":
        contents = gene_list[it].protein

        # With protein the strand direction is always positive and the positions are pulled from the object
        direc = "+"

        # Loop through each protein
        for count in range(len(contents)):

            # Pull the positions out of the protein list and turn them into strings
            start = str(contents[count].start)
            stop = str(contents[count].stop)

            # The description line is added without a number after the name
            output += "> " + gen_name(p_counter) + suffix + " " + gene_name + start + ":" + stop + ":" + direc + "\n"

            # The sequence is added
            for x in range(0, len(contents[count].sequence), 80):
                output += contents[count].sequence[x:x + 80] + '\n'

            # Increase our name counter for operon we need to add to our output
            p_counter += 1

    else:

        # The description line is added to the output string
        output += "> " + gen_name(it) + suffix + " " + gene_name + start + ":" + stop + ":" + direc + "\n"

        # The sequence is added to the output string
        for x in range(0, len(contents), 80):
            output += contents[x:x + 80] + '\n'

    return output, p_counter


def predict_helper(strand, direc, operons):
    """
    This function helps search the genome for genes and then find the RNA and protein each gene turns into

    :param strand: The DNA strand to be searched
    :param direc: Indicates if the strand is the forward (True) or reverse (False) strand
    :param operons: Determines if the strand should be searched for operons
    :return: The list of Gene objects will all attributes filled in
    """
    # gene_find is called to search the genome for genes
    # Strand is doubled so that genes that wrap around the genome can be identified
    gene_list = gene_find(strand + strand, [], direc, len(strand))

    # Each gene is transcribed, then translated
    for gene in gene_list:
        gene.rna = gene.dna.replace('T', 'U')
        translate(gene, operons)

    return gene_list


def gen_name(it):
    """
    This functions generates gene names

    :param it: A counter determining how many names have already been created
    :return: The name for the next gene
    """
    # Two strings containing the alphabet are created, they will be used to generate gene names
    l_alphabet = "abcdefghijklmnopqrstuvwxyz"
    u_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    # The genes name is generated in the format: xxxX starting with aaaA and ending with zzzZ
    # A lowercase r is appended to the gene name in the RNA file
    # A lowercase p is appended to the gene name in the Protein file
    name = l_alphabet[(it // 26 // 26 // 26) - (26 * (it // 26 // 26 // 26 // 26))] + \
           l_alphabet[(it // 26 // 26) - (26 * (it // 26 // 26 // 26))] + \
           l_alphabet[(it // 26) - (26 * (it // 26 // 26))] + \
           u_alphabet[it - (26 * (it // 26))]

    return name


main()
