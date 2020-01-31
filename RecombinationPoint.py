import queue

def read_data_as_chromosomes(file_name:str):
    f = open(file_name, 'r')
    content = f.read()
    chromosome_dict = dict()
    content = content.split("\n")
    for chromosome in content:
        chromosome_dict[chromosome[:2]] = chromosome[4:]
    return chromosome_dict

def find_mosaic_relationship():
    parent1 = input("Enter the chromosome name for the first parent chromosome: ")
    parent2 = input("Enter the chromosome name for the second parent chromosome: ")
    child = input("Enter the chromosome name for the child chromosome: ")
    return [(parent1, parent2), child]

def find_recombination_point(data, mosaic_relationship):
    parent1 = data[mosaic_relationship[0][0]]
    parent1_data = build_processed_chromosome_data(parent1)
    parent2 = data[mosaic_relationship[0][1]]
    parent2_data = build_processed_chromosome_data(parent2)
    child = data[mosaic_relationship[1]][:]
    child_data = build_processed_chromosome_data(child)
    check_queue = queue.LifoQueue()

    while len(child) > 0:

        # implement queue to track combinations from parent1 and parent2 ranges
        # and do a while loop to reduce length of child and parents for each iterative check
        first_child_char = child[0]
        parent1_ranges, parent2_ranges = compare_xdiff(first_child_char, parent1_data, parent2_data, child_data)



def compare_xdiff(char, parent1_data, parent2_data, child_data):
    parent1_xdiff = parent1_data[char+"Diff"]
    parent2_xdiff = parent2_data[char+"Diff"]
    child_xdiff = child_data[char+"Diff"]
    parent1_ranges = []
    parent2_ranges  = []
    child_counter = 0
    current_range = []

    for i in range(len(parent1_xdiff)):
        if parent1_xdiff[i] == child_xdiff[child_counter] and len(current_range) == 0:
            current_range = [1, parent1_data[char][i], i]
            child_counter += 1
        elif parent1_xdiff[i] == child_xdiff[child_counter] and len(current_range) > 0:
            current_range[2] = parent1_data[char][i + 1]
            current_range[0] = i - current_range[1]
            child_counter += 1
        else:
            if len(current_range) > 0 and current_range[0] > 1:
                parent1_ranges.append((current_range, 1))
            current_range = []
            child_counter = 0

    current_range = []
    child_counter = 0

    for i in range(len(parent2_xdiff)):
        if parent2_xdiff[i] == child_xdiff[child_counter] and len(current_range) == 0:
            current_range = [1, parent2_data[char][i], i]
            child_counter += 1
        elif parent2_xdiff[i] == child_xdiff[child_counter] and len(current_range) > 0:
            current_range[2] = parent2_data[char][i + 1]
            current_range[0] = i - current_range[1]
            child_counter += 1
        else:
            if len(current_range) > 0 and current_range[0] > 1:
                parent2_ranges.append((current_range, 2))
            current_range = []
            child_counter = 0

    parent_ranges = parent1_ranges + parent2_ranges
    parent_ranges.sort(key= lambda x: x[0][0])
    print(parent_ranges)
    return parent_ranges




def build_processed_chromosome_data(chromosome):
    chromosome_data = {"T": [], "A": [], "C":[], "G":[], "TDiff": [], "ADiff": [], "CDiff": [], "GDiff": []}
    for i in range(len(chromosome)):
        chromosome_data[chromosome[i]].append(i)
    for type in ["T", "A", "C", "G"]:
        for i in range(len(chromosome_data[type]) - 1):
            chromosome_data[type + "Diff"].append(chromosome_data[type][i + 1] - chromosome_data[type][i])
    return chromosome_data

if __name__ == "__main__":
    data = read_data_as_chromosomes("data.txt")
    mosaic_relationship = find_mosaic_relationship()
    find_recombination_point(data, mosaic_relationship)