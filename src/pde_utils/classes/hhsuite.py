import re

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
HHRESULT_HEADERS = {"query_id": re.compile("Query         ([\w\W\s\S]+)\n"),
                    "match_cols": re.compile("Match_columns ([\d]+)\n"),
                    "num_seqs": re.compile("No_of_seqs    (\d+ out of \d+)\n"),
                    "neff": re.compile("Neff          ([\d\.]+)\n"),
                    "searched_HMMs": re.compile("Searched_HMMs ([\d]+)\n"),
                    "date": re.compile("Date          ([\w\W\s\S]+)\n"),
                    "command": re.compile("Command       ([\w\W\s\S]+)\n")}

MATCHES_HEADER = (re.compile(
        " No Hit                             Prob E-value P-value  Score    "
        "SS Cols Query HMM  Template HMM\n"))

# MATCH = re.compile("\s*([\d]+)\s*([\S\-_\d\w]+)\s*(\d+\.\d+)\s*([E\w\-\d]+)\s*")

MATCH = re.compile(
         ("\s*([\d]+)\s*([\S\-_\d\w]+)\s*(\d+\.\d+)\s*([\.E\w\-\d]+)\s*"
          "([\.E\w\-\d]+)\s*(\d+\.\d+)\s*(\d+\.\d+)\s*(\d+)\s*"
          "(\d+)-(\d+)\s+(\d+)-(\d+)\s*\((\d+)\)\n"))

HHMATCH_ATTR = ["match_num", "hit_id", "probability", "expect", "p_value",
                "score", "SS", "match_cols", "query_start", "query_end",
                "hit_start", "hit_end", "hit_length"]


class HHResult:
    """Class for handling HHsuite result files.
    """
    def __init__(self, filepath):
        self.__filepath__ = filepath
        self.__initialized__ = False

        self.matches = None

        for attr_name in HHRESULT_HEADERS.keys():
            setattr(self, attr_name, None)

    def parse_result(self, file_type="hhresult"):
        """Parses HHsuite result file.
        """
        with open(self.__filepath__, "r") as filehandle:
            # Parse header of result file
            for attr_name, reg_expr in HHRESULT_HEADERS.items():
                header_line = filehandle.readline()

                header_split = re.split(reg_expr, header_line)
                if len(header_split) != 3:
                    raise HHResultFormatError("File header could not be "
                                              "identified as a HHresult file."
                                              "\nHHparser failed to recognize "
                                              f"'{attr_name}' header line in "
                                              f"{self.__filepath__}.")
                setattr(self, attr_name, header_split[1])

            filehandle.readline()
            matches_header_line = filehandle.readline()
            if re.match(MATCHES_HEADER, matches_header_line) is None:
                raise HHResultFormatError("File header could not be "
                                          "identified as a HHresult file.")

            matches = []
            while True:
                matches_line = filehandle.readline()

                if re.match(MATCH, matches_line) is None:
                    break

                matches_split = re.split(MATCH, matches_line)

                if len(matches_split) != 15:
                    raise HHResultFormatError("HHresult file matches could "
                                              "not be recognized.")

                matches_split.pop(14)
                matches_split.pop(0)

                match = HHMatch(self.query_id)
                match.load_from_list(matches_split)

                matches.append(match)

            self.matches = matches

    def check_initialization(self, caller):
        """Safe programming feature - raise an exception if a client
        tries to perform operations on uninitialized object.
        :param caller: name of the method that called this one
        """
        if not self.initialized:
            m = (f"Cannot call method '{caller}' on uninitialized "
                 "HHResult object")

            raise InitializationError(m)


class HHMatch:
    def __init__(self, query_id):
        self.query_id = query_id

        for i in range(len(HHMATCH_ATTR)):
            setattr(self, HHMATCH_ATTR[i], None)

    def load_from_list(self, match_split):
        for i in range(len(HHMATCH_ATTR)):
            setattr(self, HHMATCH_ATTR[i], match_split[i])


class InitializationError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class HHResultFormatError(Exception):
    pass
