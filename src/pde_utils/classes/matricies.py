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
        self.clusters = dict()

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

    def cluster(self, threshold):
        """
        Clusters the labels by interpreting the threshold as a lower
        bound for cluster inclusion. The matrix is temporarily cast
        to an adjacency matrix, with pairwise nodes above threshold
        being placed into the same cluster.
        :param threshold: the lower bound for cluster inclusion
        :return:
        """
        # Step 0: clear any existing clusters
        if len(self.clusters) > 0:
            self.clusters = dict()
        # Step 1: build an adjacency matrix
        adj_mat = list()
        for i in range(self.size):
            for j in range(i + 1, self.size):
                if self.matrix[i][j] >= threshold:
                    adj_mat.append({self.labels[i], self.labels[j]})
        # Step 2: build connected components from adjacency matrix
        while len(adj_mat) > 0:
            anchor = adj_mat.pop(0)
            cleanup = list()
            for i, edge in enumerate(adj_mat):
                if anchor.intersection(edge) > set():
                    anchor = anchor.union(edge)
                    cleanup.append(i)
            # Check whether this component now overlaps an existing one
            merged = False
            for key, component in self.clusters.items():
                if component.intersection(anchor) > set():
                    self.clusters[key] = component.union(anchor)
                    merged = True
                    break
            if not merged:
                self.clusters[len(self.clusters) + 1] = anchor
            for i in reversed(cleanup):
                adj_mat.pop(i)
        # Step 3: add any individual nodes that are not connected to others
        used = set()
        for key, value in self.clusters.items():
            used = used.union(value)
        missing = set(self.labels).difference(used)
        for node in missing:
            self.clusters[len(self.clusters) + 1] = {node}

    def get_clusters(self, threshold):
        """
        Calls the `cluster` method with the indicated threshold, then
        returns the resultant clusters
        :param threshold: the lower bound for cluster inclusion
        :return:
        """
        self.cluster(threshold)
        return self.clusters

    def extract_cluster_matrices(self):
        """
        Returns a list of (potentially) smaller SymmetricMatrix objects
        each of whose contents comprise one cluster from this matrix's
        data.
        Assumes this matrix has already been clustered.
        :return:
        """
        # If clustering already done, build and return cluster matrices
        if len(self.clusters) > 0:
            matrices = dict()
            for cluster_key, cluster in self.clusters.items():
                cluster = list(cluster)
                cluster_matrix = SymmetricMatrix(cluster)
                for x in range(len(cluster)):
                    label1 = cluster[x]
                    index1 = self.get_index_from_label(label1)
                    for y in range(x + 1, len(cluster)):
                        label2 = cluster[y]
                        index2 = self.get_index_from_label(label2)
                        cluster_matrix.fill_cell(x, y, self.get_cell(
                                                            index1, index2))
                cluster_matrix.fill_diagonal(100.0)
                matrices[cluster_key] = cluster_matrix
            return matrices
        # Else raise ValueError
        raise ValueError("cannot extract cluster matrices before clustering")

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
