import queue, copy

class PermutationCombination:

    def __init__(self):
        self.size = 0
        self.dict = dict()

    def add_permutation(self, keys: [], permutation, end=False):
        self.size += 1 if end else 0
        pointer = self.dict
        # TODO: change key method since None cannot equals to something
        for index, key in enumerate(keys):
            pointer = pointer[key]
        pointer[permutation] = dict()
        # print(self.dict)
    def remove_permutation(self, keys: [], permutation):
        pointer = self.dict
        for key in keys:
            pointer = pointer[key]
        # print("pre", self.dict)
        del pointer[permutation]
        # print("post", self.dict)



class RecombinationAnalysis:

    def __init__(self, file_name: str):
        data = self.read_data_as_chromosomes(file_name)
        self.parent1, self.parent2, self.child = self.map_data_to_chromosomes(data)
        self.parent1_data, self.parent2_data, self.child_data = self.build_processed_chromosome_data(self.parent1), self.build_processed_chromosome_data(self.parent2), self.build_processed_chromosome_data(self.child)
        self.analysis_result = PermutationCombination()

    def read_data_as_chromosomes(self, file_name: str):
        f = open(file_name, 'r')
        content = f.read()
        chromosome_dict = dict()
        content = content.split("\n")
        for chromosome in content:
            chromosome_dict[chromosome[:2]] = chromosome[4:]
        return chromosome_dict

    def map_data_to_chromosomes(self, data):
        parent1 = data[input("Enter the chromosome name for the first parent chromosome: ")]
        parent2 = data[input("Enter the chromosome name for the second parent chromosome: ")]
        child = data[input("Enter the chromosome name for the child chromosome: ")]
        return parent1, parent2, child

    def build_processed_chromosome_data(self, chromosome):
        chromosome_data = {"T": [], "A": [], "C": [], "G": [], "TDiff": [], "ADiff": [], "CDiff": [], "GDiff": []}
        for i in range(len(chromosome)):
            chromosome_data[chromosome[i]].append(i)
        for chromosome_type in ["T", "A", "C", "G"]:
            for i in range(len(chromosome_data[chromosome_type]) - 1):
                chromosome_data[chromosome_type + "Diff"].append(chromosome_data[chromosome_type][i + 1] - chromosome_data[chromosome_type][i])
        return chromosome_data

    def find_recombination_point(self):
        self.recursive_recombination_check(self.parent1, self.parent1_data, self.parent2, self.parent2_data, self.child, self.child_data, [])
        print(self.analysis_result.dict)
        return self.analysis_result.dict

    def recursive_recombination_check(self, parent1, parent1_data, parent2, parent2_data, child, child_data, keys: []):
        if self.analysis_result.size < 5:
            print("child", len(child))
            print("child", child[:3])
            print(keys)
            # print("parent1", len(parent1))
            # print("parent2", len(parent2))
            # print("parent2", parent2[:3])
            # print(len("GATGTCGCCTTCGAGTCCTATAGTGCGTCGCTGGTCCGCAATATTCGTCGCATCTCCATGAAGGAGGGGTGGATAGTCATCGGTCGAGATAGTAACGTAATAGTGCCTCTGACGAGCCCATGCCTAAGCGACAACCCTTGACAAAGCGCGACCTTAGGCGGACACCCCGTCCAGGCTGACTGAGGGGCAACTTTAAAAACCTACGGTACAGGCTCCTACTACGTACATACGTTTATACCCAGGTCCAGTCTTAGCGCGGTACCAGTTAATACTACACCAGTGGTCCACGCACACACTATCGAGGCGAATGTTGCCGTTAACGAAGAACCTCCCCGGAATCTGAGACCCTGTGAATGGTAAGTTAGTCCGCCAACCTCAATTATCCTTAGCTAAACACCTCGACGTAGCTCATATCACCCAGTACACTCGTAACGACTAAAAGTGTAGAGTCAAAGCAGTCCACTGTGAGGGGCAACTGGTCGTACTCTCCGGGGGAGACTTGTGTAATTAGGTGTTCGGAGGTCCCATAGTCGTCCCCTAAAGCAAGACCTCTCGCCAGGAATGTCGACGAGTTGAAGACCTGGACTTGCTACACCCGAGATAGACTAGGATGGTACTGGAGCGTAGCAGCTCATTATTGACTGCATTCAAATGCTCCGTTGGGCCTTATTTACTCGGAAAATTTACCGTGAACACTGCTAAGACTGCTCGATTCGTATTAAGGAAACATACTTTGGGTCTGAAACTATGAAACTGCACCCAGCGCCTAAAGTAGTTCAGCTAGCCTGCGCGCGGGGGGAATAGTAACTGCGTTAGCACTCTTGCAACCAGTCAACCCCTGCCCAAAGAGTGGGCTTGAGGGGGAGATCATCTTTTCAAACAATAGTTTAATTGTACTAGATCCACGACGAAGCGAAGATCGTCGATGCCTTGTTGTAACTATGACTGCTCATGCTCCGTTAAAGCCCAACGGTGTGAGGCCACCTGGCATCCCGCATGGACAATTCCTCGAGGTACAGCTCCCCATTCTCGAAATGGGCGCCGGCATGTTTACACAGGTACGTTACACTCCCTACAGTTGGCGCGTCTGGCCACAGCTAAGCCGTGGAAAAACTGTTAGCCTACGACAGTGTTGTCGCTTCTACTCGACTGGGGGCGACTGGTACCTTTAAAAATAAAGTTAGAGCAGGGGATCATCATCGATCTCGTACCTGTTGCCACAGCCTGGGCCATACCCTTAAGTGCTAGTTCCGATCTCTCAACTGGGCAGCGTGATAGACCTATTTAAATGCTAACGGCTGTGTACAAGCAGTCCAAATTGGGCCTCCTTCTTCTTGCGAGTACTAAACCCTTGACACTTTCGTAATTCACCAAACTGAACCGGGTATAATTTTCGAATAGCCGTCGCTTCTCCGGAGGACGTCCTCTTACCTCGGAGTGGGGTTTCGGCCGGACCCACGTCCATGGCCTGCTCCAAAGTTAGTTAAATCTACTTCTATTAGCAATATCCTTTATACCACAAGGGCCGAGGCTATGATGCACCACTGCAGATCACCTGAGTGATTTAAGTAATGAGGGGCAGCACGTAACGCGGCGTTATATTAGCTTCTAGAACTCCCACTCTGTGGAACGAATCTGAAATACGTCTATAGGTTGATAGCGTTCGCTAAATTTGGTCCGTACCGTACTGGGATAGTAACTAGAACATCTACCCACTATGGGAAAGGGAAGGAGCATTGGACTCCAAACTAAGGACTTTACTCACGGGATCGCGCAAGCCCAGACCAGTAGTCACCGGCTATGTCAGAGGAGACAAAAATGTCAAGCGTATTTGTTATGGTGAGACGACACAGATAGTCAAGACGGTCCAAATGTTCACCGTCGACGGACCCGATCTAATGTGAGTTACCCTTGGTGAGTGGTTCCCCGTGCGGATCGCCCTTACGCGCTGGTGAAGGGTCATCGACCCAGGAACACGGCTGAATGTAATCTTTCCCAACGCTGGCGTTCGATTTATCAGTGACAGCAGCTACAGATCCAGACGAGAAGAAGGCGATTAGAGCGTGCACGTTAGATGTTACTGACACTTCAGACGTGCACGAAGCCCCGCCTCCACGTAAGGGGTATGGGTCGAACGAGTGTCATGAGCTAACTACGAACCCTACTTTTATGGGGAAGGACCTATGTATATCCCTAAATTCAACGTCCTATGATTTGGGCACGACTATTCTCCAAACCAAAGCTTAGGCAAGGTAGCATGATACCTAACATGAACCGACTTCCCACAACGAATTTTATTATCACGTATTACTAAAGTAGAAGAGTACCTCGGCTGCGCAGATAAATTCATAACATGTCCCATTACGCCCCACGCCGGCTGCGCTAAGGAGACTTTGGATCTCCTTGCAGTGCAACATAATCGCACGCCGGTCTACACATCCGTATTTTACGGTCACCAACTCAAGATGGCTATGCACGCTATCGCAAAATCACCCCCCATTTCAAGAACGATCACAGACCCAGACATGAATCTCGTCAACATCACCAGGCCACAAGGTATTTTTCTATAAAATACTGCAGCTGAGGATAGCCGACCACGCCGGGGCAGCCAGGCCAGATTCTACTCCCATGCGACCGAACTCTATGACGGAGCGGAACCTTGTGACTTCGCAGTTTGGGGCCTATGCTCGTAGGGGTGCGAGAGCCCGAATTGGCATCTCGCGGGCCTGGGCATGCCTTTATATCCAGCCAGGTTCCGTATACGCGAACGGCGTAGTGTGCTGTCGGTGATGTCCCCGGTTCACGCGTGCGGGGATTGAGTCTGTGTGGTAGTTATTAGGCATGACACGAGCGCATATCTGAATGCCCCGGGGCCCGACCACAAATTTCTGGCTAGGCGTCCGTGCAGATGTGCGCTATTCCGTTACATAGGTACACGAGGGCTTCTGTCCATTACAAGCCCCAGAGAGTTTCCAGGCGACTCCATAATGGCATTATGGAAAGTGGGGAGCTGTTGCCGCCTCTTCTTGAAAGATTAAGGACGTCTTATGGCCGTGACATAGCAGGAAACCCCCTATCTGGCTACTACCTAGTTTCTCCTGGGTCTGAGGGGCGAGAACCGGCGGTCCTTATATTGAAACTCGGATGGTGAGGTTGCTGCTCCTACTAACATAAGCCCGATCAAAGAGTTCAGACAGGCCGTGTCGATCGACACAACAACAATCCACGTCCCGAGGGGGCGGTGCCGTACCTGGCGGCTGAATCGCTGCTGCACATTCCCCAGCTTCTAGCCGTGCGTCTCGGCTAAATTATTCCCTTTGATTCTCTCGCGCGGTCGGTTCGGACGTAGTGGAGTAGATGTTACCTAAAACGTACGAAGGCGCCGGTACTCCACGACTTGCTCTGATATCGCCCCGTCTTAACAATCAACCGCTGTTGTGAACTATAATCTACAGAAGTCACGCCACAGCGTAGGCCTGCTCATTGCGCAATGGAGGTGTGGCGTTGAGTGCTGAGGGAGGTGTTAACAAACAGGACGGGTCATTTGGATGCGAATTCCCCTCCCGCGCATGGTCCGATAAGACTAAGATGTTAACCAAGCAGGGAGATGATACTAGTTCAACCGTCAGAGACACATCTAACCAAGTAGACGAAAGGCGTGGCAGCAGCTAACCATCGAGGTAGTGAGACGTCAAAGGGTGTTAGACTTTATCAGCCAAAAGTAATATTTCTTGCGTATGAGCTTACCCTTCATGTAGGACGTAAGCGTCCAACCTTGCTACACGTAGCTGGTGCTATTTCGAGCTGCATTACCCTTCCCATTTCTTAGTAGCTGGGAAGGTTACCGTAGCCGTCACTACCTCTTACTAAGCCTAAAAGATCTTCGCGCCGTTGTAAGTTTAACCAGGAGAACCCCGCGGACTAAGTAGACGGCGGCTAGCCAATTATTTATTGCAGTGCCTTATTGGTGTCCACGTTTTCGACCTTCTAGCCGGAGTGCACTTGTTCGACTTCATCACTACCTGCAGCACTTTGAACATTGTGGACCTTCAGTCCGAGTTGTCATATTAACTGATTCGATTCTCTTACCTAGATATACTCCCAGTAACCACTTCCCCCGCCACTCTGCGGCTTATGGACAGTCAGGGCAGGACGTTCATTCGTTAATATCACTCGCCTGGTCTCATCAGTATCACAGCGAATGAGTACAGCTCCCGGGGCATCACCGACATTTATAACCCCGGTTTTTTCCAGGCCCTCTTGTAAAATTCGAACAGTGGACGCCATACCACCGGTGCAAACGGACTAATTCGTAGGACAAACATGTGGAGACAATGTTGGCTCGCAGCTAGGTGTACCCAGATAGGTAACAGCATCGAAGTATGTGACGCTCTGCGCTCTCGCCATTCCCGAATACAACATCCAGTTGTCATATCCCCCTAAAAACTCGACGGGTCAGCATAAAGGTTCTTGGGTGCGATGGGGTCCGCAGGAATTGCAAGGAGATCTTGTGCCAGGCACGCGCTGACAAGATAAACCAGGGAGATAAGGATGAGCGATTCAGCTCTAGAAGGGAAGGAACGCCGAGCGAGTAAGCCAGTGCCCCGTTTTCTGCCACGCGCCCAATTATCTTTCCAGGACAAAGAGGTCGGTGCAACTGGGGGTAGAGACTTGCCAATTCTGTATCACACATGCAGCGCAAGATTCCGGTACCTAACCCCAGGTAAGTCGATTGAGCCCCGAGGGCGGGTGTCATTGACAGTAATCGTACTGTTGGCACGCCGAACGCTCGCGCGGACTTACTGTAAATGGTGTTTGCGATACCGGGGAATGTCTTCATCACTCCTAAGACTGCGATATCGTGCCAACTGATGTCGATTTTAATATGCACTCGTTGATCTGAACGGCGGCCGCGTTCTGTTGCAGGAGTGCATGGTCGTGTAGGTGAGCACGTACTCATAATGAATGGGCAGCCATCTTGTACGAAAAGATTCACACGTACCGGCGACTGTATATCGCATATAACTTATTATCCCGGATGCGGTAGCCCCGCTGACGGCCAATAGAGGTACCACTCATTCATGTGGCTGATTCCACTT"))
            d_child = child[:]
            d_child_data = copy.deepcopy(child_data)
            d_parent1 = parent1[:]
            d_parent1_data = copy.deepcopy(parent1_data)
            d_parent2 = parent2[:]
            d_parent2_data = copy.deepcopy(parent2_data)
            first_child_char = d_child[0]
            match_dataset = self.compare_xdiff(first_child_char, d_parent1_data, d_parent2_data, d_child_data)
            print(match_dataset)
            if len(match_dataset) == 0:
                return False
            check_queue = queue.LifoQueue()
            for match_data in match_dataset:
                check_queue.put(match_data)
            while check_queue.qsize() > 0 :
                match_data = check_queue.get()
                target_parent = d_parent1 if match_data[1] == 1 else d_parent2
                target_parent_data = d_parent1_data if match_data[1] == 1 else d_parent2_data
                counter = 0
                while target_parent[match_data[0][1] + counter] == d_child[counter]:
                    parent_index = target_parent_data[target_parent[match_data[0][1] + counter]].index(match_data[0][1] + counter)
                    target_parent_data[target_parent[match_data[0][1] + counter]].pop(parent_index)
                    d_child_data[d_child[counter]].pop(0)
                    counter += 1

                recombination_point = (match_data[0][1] + counter, match_data[1])
                if match_data[1] == 1:
                    d_parent1 = d_parent1[:match_data[0][1]] + d_parent1[(match_data[0][1] + counter):]
                    d_parent1_data = self.build_processed_chromosome_data(d_parent1)
                else:
                    d_parent2 = d_parent2[:match_data[0][1]] + d_parent2[(match_data[0][1] + counter):]
                    d_parent2_data = self.build_processed_chromosome_data(d_parent2)
                d_child = d_child[counter:]
                d_child_data = self.build_processed_chromosome_data(d_child)

                if len(d_child) == 0:
                    self.analysis_result.add_permutation(keys, recombination_point, end=True)
                    return True
                else:
                    self.analysis_result.add_permutation(keys, recombination_point)
                    if self.recursive_recombination_check(d_parent1, d_parent1_data, d_parent2, d_parent2_data, d_child, d_child_data, keys + [recombination_point]):
                        return True
                    else:
                        self.analysis_result.remove_permutation(keys, recombination_point)
                        return False

    def compare_xdiff(self, char, parent1_data, parent2_data, child_data):
        parent1_xdiff = parent1_data[char+"Diff"]
        parent2_xdiff = parent2_data[char+"Diff"]
        child_xdiff = child_data[char+"Diff"]
        parent1_ranges = []
        parent2_ranges = []
        child_counter = 0
        parent_counter = 0
        current_range = []
        for i in range(len(parent1_xdiff)):
            while i + parent_counter < len(parent1_xdiff) and child_counter < len(child_xdiff) and parent1_xdiff[i + parent_counter] == child_xdiff[child_counter]:
                if len(current_range) == 0:
                    current_range.append(1)
                    current_range.append(self.parent1_data[char][i])
                    current_range.append(self.parent1_data[char][i])
                    child_counter += 1
                    parent_counter += 1
                else:
                    current_range[2] = self.parent1_data[char][i + 1]
                    current_range[0] = current_range[2] - current_range[1]
                    child_counter += 1
                    parent_counter += 1
            if len(current_range) > 0 and current_range[0] > 1:
                parent1_ranges.append((current_range, 1))
            current_range = []
            child_counter = 0
            parent_counter = 0
        print("Parent1_ranges", parent1_ranges)
        current_range = []
        child_counter = 0
        parent_counter = 0

        for i in range(len(parent2_xdiff)):
            if parent2_xdiff[i] == child_xdiff[child_counter] and len(current_range) == 0:
                current_range = [1, self.parent2_data[char][i], self.parent2_data[char][i]]
                child_counter += 1
            elif parent2_xdiff[i] == child_xdiff[child_counter] and len(current_range) > 0:
                current_range[2] = self.parent2_data[char][i + 1]
                current_range[0] = current_range[2] - current_range[1]
                child_counter += 1
            else:
                if len(current_range) > 0 and current_range[0] > 1:
                    parent2_ranges.append((current_range, 2))
                current_range = []
                child_counter = 0
        parent_ranges = parent1_ranges + parent2_ranges
        parent_ranges.sort(key=lambda x: x[0][0])
        return parent_ranges[-5:]



if __name__ == "__main__":
    RecombinationAnalysis("data.txt").find_recombination_point()