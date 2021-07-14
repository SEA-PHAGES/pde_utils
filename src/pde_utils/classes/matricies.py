import math
import numpy as np


def matrix_from_file(filepath, is_distance=False):
    """Initializes a SymmetricMatrix from a file

    :param filepath: The path to the distance matrix file
    :type filepath: pathlib.Path
    :type filepath: str
    :param is_distance: Indicates whether a cell's value is a distance metric
    :type is_distance: bool
    :return: A SymmetricMatrix object from the values in the distance matrix
    :rtype: pde_utils.classes.SymmetricMatrix
    """
    with open(filepath, "r") as mat:
        # The first line of the matrix file is the number of genes
        num_nodes = int(mat.readline())

        names = list()
        rows = list()

        # All other lines build the matrix
        for line in mat:
            # Remove trailing whitespace and split on internal whitespace
            line = line.rstrip().split()
            names.append(line[0])
            rows.append([float(x) for x in line[1:]])

    # Sanity check that we parsed the matrix properly
    if not (num_nodes == len(names) == len(rows)):
        raise Exception

    matrix = SymmetricMatrix(names, is_distance=is_distance)

    for i in range(len(names)):
        for j in range(i+1, len(names)):
            if i+1 >= len(names):
                continue

            if i == j:
                continue

            matrix.fill_cell(i, j, rows[i][j])

    return matrix


class SymmetricMatrix:
    """
    Class to facilitate storing, creating, and manipulating symmetric matrices
    (for example pairwise identity or distance matrices).
    """
    def __init__(self, labels, is_distance=False):
        """
        Constructor method for a SymmetricMatrix instance.
        """
        self.labels = labels
        self.size = len(labels)
        self.matrix = np.zeros((self.size, self.size), float)
        self.distance = is_distance

        self.mean = None
        self.median = None
        self.SD = None

    def fill_diagonal(self, value):
        """
        Fills the matrix diagonal with the input value (e.g. 1 or 100
        for identity matrices, or 0 for distance matrices).
        :param value: the value to fill the matrix diagonal with
        :return:
        """
        for i in range(self.size):
            self.fill_cell(i, i, value)

    def fill_cell(self, row, col, value):
        """
        Stores `value` in matrix[`row`][`col`] and matrix[`col`][`row`]
        :param row: the row index for locating the cell of interest
        :param col: the column index for locating the cell of interest
        :param value: the value to fill in the cell of interest
        :return:
        """
        # Check that row, col is in bounds
        if not (0 <= row < self.size and 0 <= col < self.size):
            raise ValueError(f"({row},{col}) out of bounds for "
                             f"matrix of size {self.size}")
        # If we're in bounds, always set matrix[row][col] = value
        self.matrix[row][col] = value
        # If row == col, we're done - matrix[row][col]
        # is the same cell as matrix[col][row]
        if row == col:
            return
        # Otherwise also set matrix[col][row]
        self.matrix[col][row] = value

    def get_cell(self, row, col):
        """
        Returns the value in the cell at matrix[`row`][`col`]
        :param row: the row index for locating the cell of interest
        :param col: the column index for locating the cell of interest
        :return: value
        """
        if not (0 <= row < self.size and 0 <= col < self.size):
            raise ValueError(f"({row},{col}) out of bounds for "
                             f"matrix of size {self.size}")
        return self.matrix[row][col]

    def get_row(self, index, exclude_diagonal=False):
        """
        Returns a copy of the requested matrix row.
        :param index: index of the row to be returned
        :return: copy of self.__matrix__[index]
        """
        row = self.matrix[index][:]

        if exclude_diagonal:
            row = np.delete(row, index)

        return row

    def get_index_from_label(self, label):
        """
        Search function that returns a label's index in the label list
        :param label: the label whose index should be returned
        :return: index
        """
        return self.labels.index(label)

    def get_label_from_index(self, index):
        """
        Search function that returns a label's name
        :param index: the label index to return
        :return: label
        """
        return self.labels[index]

    def get_submatrix_from_indicies(self, indicies):
        labels = []

        for index in indicies:
            labels.append(self.get_label_from_index(index))

        submatrix = SymmetricMatrix(labels)
        self.fill_submatrix(submatrix, indicies)
        return submatrix

    def get_submatrix_from_labels(self, labels):
        submatrix = SymmetricMatrix(labels)

        indicies = []
        for label in labels:
            indicies.append(self.get_index_from_label(label))

        self.fill_submatrix(submatrix, indicies)
        return submatrix

    def fill_submatrix(self, submatrix, indicies):
        for i in range(len(indicies)):
            if i == len(indicies) - 1:
                continue

            subject = indicies[i]
            for j in range(i+1, len(indicies)):
                query = indicies[j]

                value = self.matrix[subject][query]

                submatrix.matrix[i][j] = value
                submatrix.matrix[j][i] = value

    def get_centroid(self):
        """
        Interprets cells as percent identities, and returns the label
        with highest average identity to all other labels. For small
        matrices (1-2 labels), this will always be the first label.
        :return:
        """
        if 0 < self.size <= 2:
            return self.get_label_from_index(0)
        else:
            avg_pws_values = list()
            for i in range(self.size):
                row = list(self.matrix[i])
                row.pop(i)
                api = float(sum(row))/(self.size - 1)
                avg_pws_values.append(api)

            if self.distance:
                repr_value = min
            else:
                repr_value = max

            representative_index = avg_pws_values.index(
                                                    repr_value(avg_pws_values))
            return self.get_label_from_index(representative_index)

    def get_average_edge(self, label):
        """
        Calculates the average pairwise distance from the query geneid
        to every other member of the pham.
        :param geneid: the anchoring (query) gene
        :return:
        """
        index = self.get_index_from_label(label)
        values = self.get_row(index, exclude_diagonal=True)

        if len(values) < 1:
            return 0

        return float(sum(values))/len(values)

    def set_mean(self, is_distance=True):
        mean = 0
        for label in self.labels:
            mean += self.get_average_edge(label)

        mean /= self.size
        self.mean = mean
        return

        mean = 0

        if self.size <= 1:
            if is_distance:
                self.mean = 0
            else:
                self.mean = 1
            return

        edges = 0
        for i in range(self.size):
            if i == self.size - 1:
                continue

            for j in range(i+1, self.size):
                mean += self.matrix[i][j]
                edges += 1

        mean /= edges
        self.mean = mean

    def get_mean(self):
        if self.mean is None:
            self.set_mean()

        return self.mean

    def set_median(self):
        values = list()

        for i in range(self.size):
            values += list(self.get_row(i, exclude_diagonal=True))

        if len(values) <= 1:
            self.median = 0
            return

        values.sort()

        median_index = int(math.floor(len(values) / 2))
        self.median = values[median_index]

    def get_median(self):
        if self.median is None:
            self.set_median()

        return self.median

    def set_SD(self):
        std_dev = 0

        if self.size == 1:
            self.SD = std_dev
            return

        std_dev = 0
        mean = self.get_mean()

        if self.size <= 1:
            self.SD = std_dev
            return

        for i in range(self.size):
            values = self.get_row(i, exclude_diagonal=True)
            for value in values:
                std_dev += (value - mean) ** 2

        std_dev /= (self.size * (self.size - 1))
        std_dev = math.sqrt(std_dev)
        self.SD = std_dev

    def get_SD(self):
        if self.SD is None:
            self.set_SD()

        return self.SD

    def standardize(self, metric="z_score"):
        mean = self.get_mean()
        std_dev = self.get_SD()

        std_symm_matrix = SymmetricMatrix(self.labels)

        if std_dev > 0:
            std_matrix = np.zeros((self.size, self.size), float)
            for i in range(self.size):
                if i == self.size - 1:
                    continue

                for j in range(i+1, self.size):
                    value = self.matrix[i][j]
                    if metric == "z_score":
                        value = (float(value) - float(mean)) / float(std_dev)

                    std_matrix[i][j] = value
                    std_matrix[j][i] = value

            std_symm_matrix.matrix = std_matrix
        else:
            std_symm_matrix.matrix = self.matrix

        return std_symm_matrix

    def get_nearest_neighbors(self, label, threshold, is_distance=True):
        """
        Scans the query's row in the percent identity matrix to find
        all target nodes with values >= threshold. Returns the list
        in descending order (highest identity first).
        :param label: the anchoring geneid
        :type query: str
        :param threshold: the percent value threshold for inclusion
        in the return list
        :type threshold: int, float
        :return: neighbors
        :rtype: list
        """
        # Get the query row (or throw an error if the node name is invalid)
        query_index = self.get_index_from_label(label)
        query_row = self.get_row(query_index)

        # Iterate over the row - store node names where identity > threshold
        neighbor_data = list()
        for i in range(len(query_row)):
            # Skip the self-match
            if i == query_index:
                continue
            value = query_row[i]
            if is_distance:
                if value <= threshold:
                    neighbor_data.append((self.get_label_from_index(i), value))
            else:
                if value >= threshold:
                    neighbor_data.append((self.get_label_from_index(i), value))

        # Order from highest to lowest identity
        reverse = (not is_distance)
        neighbor_data.sort(key=lambda x: x[1], reverse=reverse)
        neighbors = [data[0] for data in neighbor_data]

        return neighbors

    def create_adjacency_map(self):
        adj_map = dict.fromkeys(self.labels)
        for cluster in adj_map.keys():
            adj_map[cluster] = list()

        for i in range(self.size):
            if i == self.size - 1:
                continue

            source = self.labels[i]
            for j in range(i+1, self.size):
                target = self.labels[j]
                adj_map[source].append((target, self.get_cell(i, j)))
                adj_map[target].append((source, self.get_cell(i, j)))

        return adj_map

    def __repr__(self):
        s = f"{self.size}\n"
        for i in range(self.size):
            s += (f"{self.labels[i]:<24}"
                  f"{' '.join([str(x) for x in list(self.matrix[i])])}\n")
        return s

    def __le__(self, other):
        return self.size <= other.size

    def __lt__(self, other):
        return self.size < other.size
