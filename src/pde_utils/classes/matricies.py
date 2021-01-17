import numpy as np


class SymmetricMatrix:
    """
    Class to facilitate storing, creating, and manipulating symmetric matrices
    (for example pairwise identity or distance matrices).
    """
    def __init__(self, labels):
        """
        Constructor method for a SymmetricMatrix instance.
        """
        self.labels = labels
        self.size = len(labels)
        self.matrix = np.zeros((self.size, self.size), float)

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

    def get_row(self, index):
        """
        Returns a copy of the requested matrix row.
        :param index: index of the row to be returned
        :return: copy of self.__matrix__[index]
        """
        return self.matrix[index][:]

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

                submatrix.fill_cell(i, j, self.get_cell(subject, query))

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
            avg_pws_identities = list()
            for i in range(self.size):
                row = list(self.matrix[i])
                row.pop(i)
                api = float(sum(row))/(self.size - 1)
                avg_pws_identities.append(api)
            representative_index = avg_pws_identities.index(
                                                    max(avg_pws_identities))
            return self.get_label_from_index(representative_index)

    def get_average_value(self, label):
        """
        Calculates the average pairwise distance from the query geneid
        to every other member of the pham.
        :param geneid: the anchoring (query) gene
        :return:
        """
        index = self.get_index_from_label(label)
        identities = self.get_row(index)
        identities = [identities[i] for i in range(len(identities))
                      if i != index]

        if not identities:
            return 0

        return float(sum(identities))/len(identities)

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
