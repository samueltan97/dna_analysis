import queue, copy

# PermutationCombination is the class that contains the desired answer to the possible recombination points. It also has helper methods to manage the storing and organizing of the different permutations
# for the recombination points. It basically stores a maximum of 5 most possible recombination points (which are determined as permutations with the least number of recombinations). Note how the permutations
# do not say which part of the parent chromosomes it comes from. It simply informs the user which part of the child chromosome comes from which parent chromosome.
class PermutationCombination:

    def __init__(self):
        self.size = 0
        # permutations of recombinations would be stored in self.dict. The idea is that each possible recombination point would act as a key in the dictionary which will be linked to a corresponding
        # value that will be another dictionary that will have relevant recombination points (that are possible as a result of the previous recombination point) as keys. This results in a tree-like structure
        self.dict = dict()


    def add_permutation(self, keys: [], permutation, end=False):
        # if end is true, this means this is a final recombination point that completes the child chromosome at least for this permutation. Hence, self.size has to increase by 1 to note the completion
        # of another permutation.
        self.size += 1 if end else 0
        # permutation uses the different keys stored to access the right place in self.dict to place the new recombination point (which now acts as a key to more keys if it is not the end of the child
        # chromosome).
        pointer = self.dict
        for index, key in enumerate(keys):
            pointer = pointer[key]
        pointer[permutation] = dict()

    def remove_permutation(self, keys: [], permutation):
        # TODO permutation uses the different keys stored to access the right place in self.dict to delete the recombination point that has been determined to be unable to lead to a complete child chromosome.
        pointer = self.dict
        for key in keys:
            pointer = pointer[key]
        del pointer[permutation]


# RecombinationAnalysis is the class that collects the data from the txt file and processes the information from the different chromosomes to generate the PermutationCombination which contains the results that we
# are looking for.
class RecombinationAnalysis:

    def __init__(self, file_name: str):
        # parses the text file for all chromosomes and store as data
        data = self.read_data_as_chromosomes(file_name)
        # maps data to parent1, parent2, and child by requiring input from user to determine their relationships
        self.parent1, self.parent2, self.child = self.map_data_to_chromosomes(data)
        # uses helper method (build_processed_chromosome_data) to process chromosomes to collect useful information like where the different chromosome characters appear and the gaps between them
        self.parent1_data, self.parent2_data, self.child_data = self.build_processed_chromosome_data(self.parent1), self.build_processed_chromosome_data(self.parent2), self.build_processed_chromosome_data(self.child)
        # creates a PermutationCombination object to store possible recombination points
        self.analysis_result = PermutationCombination()

    # read_data_as_chromosomes opens the text file and builds a dictionary that contains all of the chromosomes with their names as keys
    def read_data_as_chromosomes(self, file_name: str):
        f = open(file_name, 'r')
        content = f.read()
        chromosome_dict = dict()
        content = content.split("\n")
        for chromosome in content:
            chromosome_dict[chromosome[:2]] = chromosome[4:]
        return chromosome_dict

    # requires input from user to determine which chromosomes are the parent chromosomes and the child chromosome of interest
    def map_data_to_chromosomes(self, data):
        parent1 = data[input("Enter the chromosome name for the first parent chromosome: ")]
        parent2 = data[input("Enter the chromosome name for the second parent chromosome: ")]
        child = data[input("Enter the chromosome name for the child chromosome: ")]
        return parent1, parent2, child

    # creates a dictionary (for target chromosome) that contains the positions of each T, A, C, and G, as well as the number of characters which separate each subsequent similar characters (TDiff, ADiff, CDiff, and GDiff)
    def build_processed_chromosome_data(self, chromosome):
        chromosome_data = {"T": [], "A": [], "C": [], "G": [], "TDiff": [], "ADiff": [], "CDiff": [], "GDiff": []}
        for i in range(len(chromosome)):
            chromosome_data[chromosome[i]].append(i)
        for chromosome_type in ["T", "A", "C", "G"]:
            for i in range(len(chromosome_data[chromosome_type]) - 1):
                chromosome_data[chromosome_type + "Diff"].append(chromosome_data[chromosome_type][i + 1] - chromosome_data[chromosome_type][i])
        return chromosome_data

    # critical function that needs to be called which will return the self.dict in the PermutationCombination object. It initiates a recursive function that would modify the PermutationCombination object to return
    # the possible permutations of recombination points
    def find_recombination_point(self):
        self.recursive_recombination_check(self.parent1, self.parent1_data, self.parent2, self.parent2_data, self.child, self.child_data, [])
        return self.analysis_result.dict

    # this function takes in the parent chromosomes and child chromosome, as well as their processed data along with a set of keys (the recombination permutation that it belongs to: this facilitates the modification of the right node in
    # the PermutationCombination object). It primarily stores the possible permutations in a stack before calling itself recursively to check the possible permutations and appending these possible permutattions to the PermutationCombination
    # object while removing all the failed combination checks from the PermutationCombination object. It compares the pattern of XDiff of both the parent chromosome and the child chromosome where X is the first character of the child (or
    # rather the remaining child passed to the recursive function).
    def recursive_recombination_check(self, parent1, parent1_data, parent2, parent2_data, child, child_data, keys: []):
        # TODO if the PermutationCombination object aready contains 5 optimal permutations, terminate this recursive search
        if self.analysis_result.size < 5:
            # print(len("GATGTCGCCTTCGAGTCCTATAGTGCGTCGCTGGTCCGCAATATTCGTCGCATCTCCATGAAGGAGGGGTGGATAGTCATCGGTCGAGATAGTAACGTAATAGTGCCTCTGACGAGCCCATGCCTAAGCGACAACCCTTGACAAAGCGCGACCTTAGGCGGACACCCCGTCCAGGCTGACTGAGGGGCAACTTTAAAAACCTACGGTACAGGCTCCTACTACGTACATACGTTTATACCCAGGTCCAGTCTTAGCGCGGTACCAGTTAATACTACACCAGTGGTCCACGCACACACTATCGAGGCGAATGTTGCCGTTAACGAAGAACCTCCCCGGAATCTGAGACCCTGTGAATGGTAAGTTAGTCCGCCAACCTCAATTATCCTTAGCTAAACACCTCGACGTAGCTCATATCACCCAGTACACTCGTAACGACTAAAAGTGTAGAGTCAAAGCAGTCCACTGTGAGGGGCAACTGGTCGTACTCTCCGGGGGAGACTTGTGTAATTAGGTGTTCGGAGGTCCCATAGTCGTCCCCTAAAGCAAGACCTCTCGCCAGGAATGTCGACGAGTTGAAGACCTGGACTTGCTACACCCGAGATAGACTAGGATGGTACTGGAGCGTAGCAGCTCATTATTGACTGCATTCAAATGCTCCGTTGGGCCTTATTTACTCGGAAAATTTACCGTGAACACTGCTAAGACTGCTCGATTCGTATTAAGGAAACATACTTTGGGTCTGAAACTATGAAACTGCACCCAGCGCCTAAAGTAGTTCAGCTAGCCTGCGCGCGGGGGGAATAGTAACTGCGTTAGCACTCTTGCAACCAGTCAACCCCTGCCCAAAGAGTGGGCTTGAGGGGGAGATCATCTTTTCAAACAATAGTTTAATTGTACTAGATCCACGACGAAGCGAAGATCGTCGATGCCTTGTTGTAACTATGACTGCTCATGCTCCGTTAAAGCCCAACGGTGTGAGGCCACCTGGCATCCCGCATGGACAATTCCTCGAGGTACAGCTCCCCATTCTCGAAATGGGCGCCGGCATGTTTACACAGGTACGTTACACTCCCTACAGTTGGCGCGTCTGGCCACAGCTAAGCCGTGGAAAAACTGTTAGCCTACGACAGTGTTGTCGCTTCTACTCGACTGGGGGCGACTGGTACCTTTAAAAATAAAGTTAGAGCAGGGGATCATCATCGATCTCGTACCTGTTGCCACAGCCTGGGCCATACCCTTAAGTGCTAGTTCCGATCTCTCAACTGGGCAGCGTGATAGACCTATTTAAATGCTAACGGCTGTGTACAAGCAGTCCAAATTGGGCCTCCTTCTTCTTGCGAGTACTAAACCCTTGACACTTTCGTAATTCACCAAACTGAACCGGGTATAATTTTCGAATAGCCGTCGCTTCTCCGGAGGACGTCCTCTTACCTCGGAGTGGGGTTTCGGCCGGACCCACGTCCATGGCCTGCTCCAAAGTTAGTTAAATCTACTTCTATTAGCAATATCCTTTATACCACAAGGGCCGAGGCTATGATGCACCACTGCAGATCACCTGAGTGATTTAAGTAATGAGGGGCAGCACGTAACGCGGCGTTATATTAGCTTCTAGAACTCCCACTCTGTGGAACGAATCTGAAATACGTCTATAGGTTGATAGCGTTCGCTAAATTTGGTCCGTACCGTACTGGGATAGTAACTAGAACATCTACCCACTATGGGAAAGGGAAGGAGCATTGGACTCCAAACTAAGGACTTTACTCACGGGATCGCGCAAGCCCAGACCAGTAGTCACCGGCTATGTCAGAGGAGACAAAAATGTCAAGCGTATTTGTTATGGTGAGACGACACAGATAGTCAAGACGGTCCAAATGTTCACCGTCGACGGACCCGATCTAATGTGAGTTACCCTTGGTGAGTGGTTCCCCGTGCGGATCGCCCTTACGCGCTGGTGAAGGGTCATCGACCCAGGAACACGGCTGAATGTAATCTTTCCCAACGCTGGCGTTCGATTTATCAGTGACAGCAGCTACAGATCCAGACGAGAAGAAGGCGATTAGAGCGTGCACGTTAGATGTTACTGACACTTCAGACGTGCACGAAGCCCCGCCTCCACGTAAGGGGTATGGGTCGAACGAGTGTCATGAGCTAACTACGAACCCTACTTTTATGGGGAAGGACCTATGTATATCCCTAAATTCAACGTCCTATGATTTGGGCACGACTATTCTCCAAACCAAAGCTTAGGCAAGGTAGCATGATACCTAACATGAACCGACTTCCCACAACGAATTTTATTATCACGTATTACTAAAGTAGAAGAGTACCTCGGCTGCGCAGATAAATTCATAACATGTCCCATTACGCCCCACGCCGGCTGCGCTAAGGAGACTTTGGATCTCCTTGCAGTGCAACATAATCGCACGCCGGTCTACACATCCGTATTTTACGGTCACCAACTCAAGATGGCTATGCACGCTATCGCAAAATCACCCCCCATTTCAAGAACGATCACAGACCCAGACATGAATCTCGTCAACATCACCAGGCCACAAGGTATTTTTCTATAAAATACTGCAGCTGAGGATAGCCGACCACGCCGGGGCAGCCAGGCCAGATTCTACTCCCATGCGACCGAACTCTATGACGGAGCGGAACCTTGTGACTTCGCAGTTTGGGGCCTATGCTCGTAGGGGTGCGAGAGCCCGAATTGGCATCTCGCGGGCCTGGGCATGCCTTTATATCCAGCCAGGTTCCGTATACGCGAACGGCGTAGTGTGCTGTCGGTGATGTCCCCGGTTCACGCGTGCGGGGATTGAGTCTGTGTGGTAGTTATTAGGCATGACACGAGCGCATATCTGAATGCCCCGGGGCCCGACCACAAATTTCTGGCTAGGCGTCCGTGCAGATGTGCGCTATTCCGTTACATAGGTACACGAGGGCTTCTGTCCATTACAAGCCCCAGAGAGTTTCCAGGCGACTCCATAATGGCATTATGGAAAGTGGGGAGCTGTTGCCGCCTCTTCTTGAAAGATTAAGGACGTCTTATGGCCGTGACATAGCAGGAAACCCCCTATCTGGCTACTACCTAGTTTCTCCTGGGTCTGAGGGGCGAGAACCGGCGGTCCTTATATTGAAACTCGGATGGTGAGGTTGCTGCTCCTACTAACATAAGCCCGATCAAAGAGTTCAGACAGGCCGTGTCGATCGACACAACAACAATCCACGTCCCGAGGGGGCGGTGCCGTACCTGGCGGCTGAATCGCTGCTGCACATTCCCCAGCTTCTAGCCGTGCGTCTCGGCTAAATTATTCCCTTTGATTCTCTCGCGCGGTCGGTTCGGACGTAGTGGAGTAGATGTTACCTAAAACGTACGAAGGCGCCGGTACTCCACGACTTGCTCTGATATCGCCCCGTCTTAACAATCAACCGCTGTTGTGAACTATAATCTACAGAAGTCACGCCACAGCGTAGGCCTGCTCATTGCGCAATGGAGGTGTGGCGTTGAGTGCTGAGGGAGGTGTTAACAAACAGGACGGGTCATTTGGATGCGAATTCCCCTCCCGCGCATGGTCCGATAAGACTAAGATGTTAACCAAGCAGGGAGATGATACTAGTTCAACCGTCAGAGACACATCTAACCAAGTAGACGAAAGGCGTGGCAGCAGCTAACCATCGAGGTAGTGAGACGTCAAAGGGTGTTAGACTTTATCAGCCAAAAGTAATATTTCTTGCGTATGAGCTTACCCTTCATGTAGGACGTAAGCGTCCAACCTTGCTACACGTAGCTGGTGCTATTTCGAGCTGCATTACCCTTCCCATTTCTTAGTAGCTGGGAAGGTTACCGTAGCCGTCACTACCTCTTACTAAGCCTAAAAGATCTTCGCGCCGTTGTAAGTTTAACCAGGAGAACCCCGCGGACTAAGTAGACGGCGGCTAGCCAATTATTTATTGCAGTGCCTTATTGGTGTCCACGTTTTCGACCTTCTAGCCGGAGTGCACTTGTTCGACTTCATCACTACCTGCAGCACTTTGAACATTGTGGACCTTCAGTCCGAGTTGTCATATTAACTGATTCGATTCTCTTACCTAGATATACTCCCAGTAACCACTTCCCCCGCCACTCTGCGGCTTATGGACAGTCAGGGCAGGACGTTCATTCGTTAATATCACTCGCCTGGTCTCATCAGTATCACAGCGAATGAGTACAGCTCCCGGGGCATCACCGACATTTATAACCCCGGTTTTTTCCAGGCCCTCTTGTAAAATTCGAACAGTGGACGCCATACCACCGGTGCAAACGGACTAATTCGTAGGACAAACATGTGGAGACAATGTTGGCTCGCAGCTAGGTGTACCCAGATAGGTAACAGCATCGAAGTATGTGACGCTCTGCGCTCTCGCCATTCCCGAATACAACATCCAGTTGTCATATCCCCCTAAAAACTCGACGGGTCAGCATAAAGGTTCTTGGGTGCGATGGGGTCCGCAGGAATTGCAAGGAGATCTTGTGCCAGGCACGCGCTGACAAGATAAACCAGGGAGATAAGGATGAGCGATTCAGCTCTAGAAGGGAAGGAACGCCGAGCGAGTAAGCCAGTGCCCCGTTTTCTGCCACGCGCCCAATTATCTTTCCAGGACAAAGAGGTCGGTGCAACTGGGGGTAGAGACTTGCCAATTCTGTATCACACATGCAGCGCAAGATTCCGGTACCTAACCCCAGGTAAGTCGATTGAGCCCCGAGGGCGGGTGTCATTGACAGTAATCGTACTGTTGGCACGCCGAACGCTCGCGCGGACTTACTGTAAATGGTGTTTGCGATACCGGGGAATGTCTTCATCACTCCTAAGACTGCGATATCGTGCCAACTGATGTCGATTTTAATATGCACTCGTTGATCTGAACGGCGGCCGCGTTCTGTTGCAGGAGTGCATGGTCGTGTAGGTGAGCACGTACTCATAATGAATGGGCAGCCATCTTGTACGAAAAGATTCACACGTACCGGCGACTGTATATCGCATATAACTTATTATCCCGGATGCGGTAGCCCCGCTGACGGCCAATAGAGGTACCACTCATTCATGTGGCTGATTCCACTT"))
            # print(len("TAACGGGACAGAATTCGTTACATTATAGGCCGGCTGATATCCTTAGACAATCTTGCCGTACTCCAGGGCGGTGGTGGCCCCACGGCATAGGGTGGATCGGACCTAATTCTAAATCGTTAAAAATAGCCGGGCTGTGCTGGCGGCACTAGATCCGATTCATAGGTCGTTAGCCGGCGACCAGCACCTTCCGTATGCATGAAAGCGTCAGGAATGAAACCGGCAGGCTTGTAAGCCCAGTACAAGTCAGCCAGAACCAAGTAGCGGGGCCTAATGCCACGGAATCCAATTCAGCTAGCTCCTGCACTTTCCGCGCCCGATGTCTGTTTTATCGCGGTGAGAAGTGAGGCGTGTCATTACAAGGGATTTAATGTTGCTGGCCTTGAGGCCCCCCGCCGTTAGCAACTTAAAGCGTTAGGTAGTCACGTTCGACAGGCTTGAGAAATTTGATAACTCGTAGCATCTTGCGTTCCCGGGCGCAATCACGCTCACTGTAACATAATGGGTCTAATGTGACAAGCTCATCAACGCTCCCCATTTATATGCAAAATGGGACGACAAGTGCCACTACGTCCCGTTGGCGTCCGAGCAATTTTATCGATGTCTCTCTCGGAAACTAAGACTACAGGTCGAGAAGACAGGAATGTGTCCTTCGACTGTCGTTCCTAGATACGCCTATTCTGGCAGGCCCCGCTAGGTCGTTAGACAATTGCAGTTCGACATATCGGCTACACGTAGCAACGTCTTGGGATCGAACGCCGCGAAGCCTCCCCTCGAGTCAGCCTGATGTTCATAAGTAACCTAATAGGATCCTACACCGCTGGAGCGGTAAATTTTAAGCTTACAGGGCCGGATGCCTAGAGGATCACTATAGGCGAGTGGATGAGTATCAATGTTGCTTAATGATTGGCCGGAAAACTTGTAATGCAATTCTATCCTGCGAAGCAAGGGCGCCTAAACAAACGTCCAGTGCTAGAAATTGTTTCTCCTGGCTATGCATAACCAGGCTAGGCAAACAATGCTTTCTCGCTTAGTACATTACGAGGTTCTTGATCCTTGTGCCGGATGCGCCAAATTCGTCGTACACGTGGATGCACCAGACGCTGTTCTGACGTGCCCAGAAGAGATTAACACCTAGCCCGAATTTGAGACACTGGACAGAGTGCCGAAGCTAATTTGTCTGGCCAACGCAGACACCTGGCGGATTTTGTCAGCTAACACCTGCTATGCTAGATTACAATGACGGACACTTGGAGTGAGGTGAACTTACTCATGAGTGCTTACACCACCAATAAGTACGCCATGGTCATGAACTGCGTGGAGGCGAGTAAGTGTTCTTTAATGGTGCCCGCTAAGCGTGAAAGTGCCTCGTCGTCTAGTCTGTAGGGACCTGGGGTTTGTATTAGAGCATATTCCAGGAGACTCGCTTTACTACTTTTCTTTACGAGTGCTGTTCCCGTCATCCCAATTCCCAGAGGTCGCCTCGAGTCTCCATCATGGGTCCGTTTCCAGTCGCAGTTACCTGCCCGAAGGGCCAGTCTTAATGTCTCCGTAATCGCTCTTGTCCCAACTGGCCTCATCACTGATTCGCCGGGAAGCAACGCTCACACCCTGGTGCACGTACGAATTCCCGCGATAGGCATGACTTGCGGCAAGTTAGAGCATAATATTTAGTACCACAGTCGGTGGTACTGCAATAACCTTTTTTATAAGGACTAGTTCACCATAGCTATCGTCCAGCCACGTTCCGGGTGAATTCGAGGGCCCTCAAGCGCAGCTGTCGACTAGTTTATATAAACATTCTTTTACAGCTAGATGGCTCGTTCTGGATGGGGCTCAACCAACAGTATAGTCAAAAACCAGACCCATACTTAGTGCTGAACACAGATGTAAAACGAAGTGGGCGTGGAAATCGCAAAGCCGGATACGAAGACGTTGGGCCCCAGTCACAAAGTTAAACCCCTAGACGTCAAAGAGCGCCTTATGAGCTGGAAACCTTGATAGAGACTAGTAATGACTCGTACGGTGTATTCCAGGTCTAGCGAGTTGCTGGGCTTTCTCTTAAGTAATTTTGGCTACGACCGACTAGTTGATTTCCTTGCCAATCGTGGTGGTACTTGAGAACCCGTCGATAGATCATTCGACTTATACTGTAGCGCAGGGAGCTACGTAGGCTCTTGATCTCTCTACCTGACGTCCAAAATGGGCATGATGCGAAGTCCGAGGAATCATTAGTAAAACTTAGTGGTGTCACTCCCCCGGGTTAGGCTACCCCGACGTCAGTCTAGAGATGTGTCCAGCGACATTCGAAGGAGGTGTAACGTGCATTATTTTGGCGAGGCAGCAGAACATTCACCCTTACTTGATTACTTCATCCCAGGATGTCCTCCCACACTAAGACATCCGATATACTAATCGGTTCCGTGGCTATGTCGCAGAGTTTTTACCCCCGTTATGACAGTGAACCTAGTTTACAATCTGTGCGAAGTGTGATGACTCGAAGAAACCAGCAGAGGTCGGCGGAGTGGCTCGCTCATAAAGGGAAAATACTTTCGTATACAAATTACAATTTGACGATCGATTAAGCCGGTTCTACCGTACTCCACCGTATCGGGGAAATGAACAGGATATACCTGCTTAACAACAGCTTTTGCGTAGTTGGAACAAACACAGCACCTAGCCTTCAGATTGGGCGTCTAACAGCAGGATCAGTCTCCCCCATTGGATCATCGTCGCAACGCTGGTCGCCGGCAACAAAACGACATCGCAATTTGGCCGGGACGGATATTTCGTTATCAATAACGGATGTAGAAAAGACTTGGGAATACCGTTTTTGAATTGATTGTCTATCCGCTGCGTCCGACATCATCGTAAAACCAACCACCGTAAAGACATACTACTGTTGGTGCACTCGCATAGCCCCGAGGATCCTCTCTGAGCCCGCGACTTGCCTTGCATACACAAAAGAGTGGGCTGTATATTGAACCTATGGGGTAATCAGAGGACGGCCACGTAATCTGTTGATAGAAAAACTCGTTCTTTTGGCCACGAAGGGGCCATGAGGCCGCTATATCCTGCGCCTACTTTGCGGTCCAGGTGGTCCGAGCCTAAATGTGATCAGCCGATCAGACTTATTCTAAGCTATGATACTCGGTTAATATTCATTGAATATAACCTGAGATTGCGACCTCGGAGATGGTCTGCCTGCGTTCGGATGATGACAGAAGAGGGT"))

            # deep copy every single parent and child chromosomes and processed data since we will be doing modifications to them and passing the modified information to the recursive function to use them to search for the subsequent recombination points
            # in this permutation
            d_child = child[:]
            d_child_data = copy.deepcopy(child_data)
            d_parent1 = parent1[:]
            d_parent1_data = copy.deepcopy(parent1_data)
            d_parent2 = parent2[:]
            d_parent2_data = copy.deepcopy(parent2_data)
            # first_child_char is the first character of the child passed to the function
            first_child_char = d_child[0]
            # match_dataset is a list of the top 5 XDiff pattern match between the child chromosome and parent chromosomes
            match_dataset = self.compare_xdiff(first_child_char, d_parent1_data, d_parent2_data, d_child_data)
            # If there is no match at all, False will be returned to terminate the search for possible recombination points under the premise of this permutation
            if len(match_dataset) == 0:
                return False
            # A stack called check_queue is created to store the 5 most likely pattern match that would then be iteratively checked to determine if they would eventually result in a permutation that completes the child_chromosome
            check_queue = queue.LifoQueue()
            for match_data in match_dataset:
                check_queue.put(match_data)
            while check_queue.qsize() > 0 :
                # use the match_data by popping the stack. The match_data contains the following information: the number of matches of XDiff with the child_chromosome, the start and end index of parent_chromosome that matches
                # by virtue of XDiff (meaning the match could be longer just that the subsequent characters are not of X-type), and which parent chromosome are we working with now
                match_data = check_queue.get()
                # target_parent is set to point at the right parent_chromosome copy
                target_parent = d_parent1 if match_data[1] == 1 else d_parent2
                # target_parent_data is set to point at the right parent_chromosome_data copy
                target_parent_data = d_parent1_data if match_data[1] == 1 else d_parent2_data
                # counter is supposed to track the actual number of matches between target parent and child chromosome
                counter = 0
                #use a while loop to compare parent_chromosome segment with child_chromosome. If there is a match, remove that character from the processed data of both parent and child chromosome
                while match_data[0][1] + counter < len(target_parent) and counter < len(d_child) and target_parent[match_data[0][1] + counter] == d_child[counter]:
                    parent_index = target_parent_data[target_parent[match_data[0][1] + counter]].index(match_data[0][1] + counter)
                    target_parent_data[target_parent[match_data[0][1] + counter]].pop(parent_index)
                    d_child_data[d_child[counter]].pop(0)
                    counter += 1
                # TODO by using the final value of the counter, we can determine where the recombination was done in the child_chromosome as well as from which parent chromosome it came from
                recombination_point = (match_data[0][1] + counter - 1, match_data[1])
                # if the chromosome segment came from parent 1 or 2, remove the segment from the parent chromosome copy and rebuild a set of processed data for the parent chromosome copy
                if match_data[1] == 1:
                    d_parent1 = d_parent1[:match_data[0][1]] + d_parent1[(match_data[0][1] + counter):]
                    d_parent1_data = self.build_processed_chromosome_data(d_parent1)
                else:
                    d_parent2 = d_parent2[:match_data[0][1]] + d_parent2[(match_data[0][1] + counter):]
                    d_parent2_data = self.build_processed_chromosome_data(d_parent2)
                # remove matched chromosome segment from the child chromosome copy
                d_child = d_child[counter:]
                d_child_data = self.build_processed_chromosome_data(d_child)

                # if child_chromosome copy is now empty, it means that we have managed to complete a permutation for the recombination points and we will hence add that permutation to the
                # PermutationCombination object with end set as True
                if len(d_child) == 0:
                    self.analysis_result.add_permutation(keys, recombination_point, end=True)
                    return True
                # else we will just add the permutation while recursively calling this same function again with the modified copies of the chromosomes and processed data
                else:
                    self.analysis_result.add_permutation(keys, recombination_point)
                    if self.recursive_recombination_check(d_parent1, d_parent1_data, d_parent2, d_parent2_data, d_child, d_child_data, keys + [recombination_point]):
                        return True
                    # if there is no match or if we already have 5 permutations, False will be returned and we will remove the permutation from the PermutationCombination object
                    # (think of it like backtracking without accidentally clearing up the correct permutation data from the tree)
                    else:
                        self.analysis_result.remove_permutation(keys, recombination_point)
                        return False
        else:
            return False

    # this function simply compares parent1_data and parent2_data against the child_data with respect to the first character in child_data. It tries to find the 5 longest matches in XDiff (the gap between the As, Ts, Cs, or Gs depending
    # on what the first character of the child is that is denoted by char) in the different possible matches for parent1 and parent2 with child.
    def compare_xdiff(self, char, parent1_data, parent2_data, child_data):
        parent1_xdiff = parent1_data[char+"Diff"]
        parent2_xdiff = parent2_data[char+"Diff"]
        child_xdiff = child_data[char+"Diff"]
        parent1_ranges = []
        parent2_ranges = []
        child_counter = 0
        parent_counter = 0
        # current_range will be presented in this form [number of matches in xdiff, starting index of match in parent, last index of match in parent that ends with char]
        current_range = []
        # iterates through parent1_xdiff and child_xdiff in a sliding window fashion to find matches
        for i in range(len(parent1_xdiff)):
            # a while loop is used here because we want to retain the value of i to make sure we use every xdiff in parent1 as a possible starting point. child_counter and parent_counter will be used to manage the "for-loop"
            # within the actual for-loop
            while i + parent_counter < len(parent1_xdiff) and child_counter < len(child_xdiff) and parent1_xdiff[i + parent_counter] == child_xdiff[child_counter]:
                # if this is the first match between parent1_xdiff and child_xdiff, the following information will be appended into current_range
                if len(current_range) == 0:
                    current_range.append(1)
                    current_range.append(parent1_data[char][i])
                    current_range.append(parent1_data[char][i])
                    child_counter += 1
                    parent_counter += 1
                # if this is not the first match, the values in current_range will be modified to reflect the new last index of match as well as the number of matches
                else:
                    current_range[2] = parent1_data[char][i + parent_counter + 1]
                    current_range[0] = current_range[2] - current_range[1] + 1
                    child_counter += 1
                    parent_counter += 1
            # if there are no more matches, we append current range as well as 1 which stands for parent 1 into parent1_ranges before resetting current_range, and the counters
            if len(current_range) > 0 and current_range[0] > 1:
                parent1_ranges.append((current_range, 1))
            current_range = []
            child_counter = 0
            parent_counter = 0

        # the above is repeated for parent2_xdiff
        for i in range(len(parent2_xdiff)):
            while i + parent_counter < len(parent2_xdiff) and child_counter < len(child_xdiff) and parent2_xdiff[i + parent_counter] == child_xdiff[child_counter]:
                if len(current_range) == 0:
                    current_range.append(1)
                    current_range.append(parent2_data[char][i])
                    current_range.append(parent2_data[char][i])
                    child_counter += 1
                    parent_counter += 1
                else:
                    current_range[2] = parent2_data[char][i + parent_counter + 1]
                    current_range[0] = current_range[2] - current_range[1] + 1
                    child_counter += 1
                    parent_counter += 1
            if len(current_range) > 0 and current_range[0] > 1:
                parent2_ranges.append((current_range, 2))
            current_range = []
            child_counter = 0
            parent_counter = 0

        # the parent ranges are then put together and sorted by maximum number of matches. The five with the most number of matches are returned to the recursive function to be used as match_data
        parent_ranges = parent1_ranges + parent2_ranges
        parent_ranges.sort(key=lambda x: x[0][0])
        return parent_ranges[-5:]



if __name__ == "__main__":
    print(RecombinationAnalysis("data.txt").find_recombination_point())