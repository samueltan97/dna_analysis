import queue, copy

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
    recursive_recombination_check(parent1, parent1_data, parent2, parent2_data, child, child_data)

def recursive_recombination_check(parent1, parent1_data, parent2, parent2_data, child, child_data):
    d_child = child[:]
    d_child_data = copy.deepcopy(child_data)
    d_parent1 = parent1[:]
    d_parent1_data = copy.deepcopy(parent1_data)
    d_parent2 = parent2[:]
    d_parent2_data = copy.deepcopy(parent2_data)
    first_child_char = d_child[0]
    parent_ranges = compare_xdiff(first_child_char, d_parent1_data, d_parent2_data, d_child_data)
    check_queue = queue.LifoQueue()
    possible_recombination_points = []
    for parent_range in parent_ranges:
        check_queue.put(parent_range)
    while check_queue.qsize() > 0:
        parent_range = check_queue.get()
        target_parent = (d_parent1, d_parent1_data) if parent_range[1] == 1 else (d_parent2, d_parent2_data)
        counter = 0
        while (target_parent[0][parent_range[0][1] + counter] == d_child[counter]):
            parent_index = target_parent[1][target_parent[0][parent_range[0][1] + counter]].index(parent_range[0][1] + counter)
            target_parent[1][target_parent[0][parent_range[0][2] + 1 + counter]].pop(parent_index)
            d_child_data[d_child[counter]].pop(0)
            counter += 1
        # figure out how to append (basically we are only tracking which indexes of child belong to which parent)
        for chromosome_diff in ["TDiff", "ADiff", "CDiff", "GDiff"]:
            d_child_data[chromosome_diff] = []
            target_parent[1][chromosome_diff] = []

        for chromosome_type in ["T", "A", "C", "G"]:
            for i in range(len(target_parent[1][chromosome_type]) - 1):
                target_parent[1][chromosome_type + "Diff"].append(
                    target_parent[1][chromosome_type][i + 1] - target_parent[1][chromosome_type][i])
            for i in range(len(d_child_data[chromosome_type]) - 1):
                d_child_data[chromosome_type + "Diff"].append(
                    d_child_data[chromosome_type][i + 1] - d_child_data[chromosome_type][i])

        possible_recombination_points.append((parent_range[0][2] + counter, parent_range[1]))
        if parent_range[1] == 1:
            d_parent1 = d_parent1[:parent_range[0][1]] + d_parent1[(parent_range[0][1] + counter + 1):]
        else:
            d_parent2 = d_parent2[:parent_range[0][1]] + d_parent2[(parent_range[0][1] + counter + 1):]
        d_child = d_child[counter+1:]

        possible_recombination_points += recursive_recombination_check(d_parent1, d_parent1_data, d_parent2, d_parent2_data, d_child, d_child_data)
        return possible_recombination_points





def compare_xdiff(char, parent1_data, parent2_data, child_data):
    parent1_xdiff = parent1_data[char+"Diff"]
    parent2_xdiff = parent2_data[char+"Diff"]
    child_xdiff = child_data[char+"Diff"]
    parent1_ranges = []
    parent2_ranges = []
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
    parent_ranges.sort(key=lambda x: x[0][0])
    print(parent_ranges)
    return parent_ranges




def build_processed_chromosome_data(chromosome):
    chromosome_data = {"T": [], "A": [], "C": [], "G": [], "TDiff": [], "ADiff": [], "CDiff": [], "GDiff": []}
    for i in range(len(chromosome)):
        chromosome_data[chromosome[i]].append(i)
    for chromosome_type in ["T", "A", "C", "G"]:
        for i in range(len(chromosome_data[chromosome_type]) - 1):
            chromosome_data[chromosome_type + "Diff"].append(chromosome_data[chromosome_type][i + 1] - chromosome_data[chromosome_type][i])
    return chromosome_data

if __name__ == "__main__":
    data = read_data_as_chromosomes("data.txt")
    mosaic_relationship = find_mosaic_relationship()
    find_recombination_point(data, mosaic_relationship)