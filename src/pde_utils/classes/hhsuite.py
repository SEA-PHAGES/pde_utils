import re

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------

# Initial HHresult header identification regular expressions
# ==========================================================
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
# ==========================================================

# HHresult data table parsing regular expressions
# ===============================================
# Parses one HHresult hit line start into two capturing groups: No and Hit
TABLE_INDEX = re.compile(
         ("\s*([\d]+)"  # 'No'/'match_num'
          "\s*([\S\s\-_\d\w]+)$"))  # 'HitID'/'hit_id'

# Parsers one HHresult hit line tail into eleven capturing groups
TABLE_DATA = re.compile(
         ("\s*(\d+\.\d+)"  # 'Prob'/'probability'
          "\s*([\.E\w\-\d]+)"  # 'E-value'/'expect'
          "\s*([\.E\w\-\d]+)"  # 'P-value'/'p_value'
          "\s*(\d+\.\d+)"  # 'Score'/'score'
          "\s*(\d+\.\d+)"  # 'SS' secondary structure
          "\s*(\d+)"  # 'Cols'/'match_cols'
          "\s*(\d+)-(\d+)"  # 'Query HMM'/('query_start')-('query_end')
          # 'Template HMM'/('hit_start')-('hit_end') ('hit_length)
          "\s+(\d+)-(\d+)\s*\((\d+)\)\n"))
# ===============================================

# HHresult data body parsing regular expressions
# ==============================================
# Parses one HHresult body hit line start into one capturing group
BODY_INDEX = re.compile("No (\d+)")  # 'No'/'match_num'
BODY_NAME = re.compile(">([\S\s\-_\d\w]+)\n")

# Parses one HHresult body hit line start into eight capturing groups
BODY_DATA = re.compile(
                ("Probab=([\d\.]+)"  # 'probability'
                 "\s*E-value=([\d\.\w\-]+)"  # 'expect'
                 "\s*Score=([\d\.]+)"  # 'score'
                 "\s*Aligned_cols=(\d+)"  # 'match_cols"
                 "\s*Identities=([\d\.]+)%"  # 'pid'
                 "\s*Similarity=([\d\.]+)"  # 'similarity'
                 "\s*Sum_probs=([\d\.]+)"  # 'sum_probs'
                 "\s*Template_Neff=([\d\.]+)"))  # 'template_Neff'

# Parses one HHresult body alignment line into three capturing groups
BODY_ALIGNMENT_SEQ = re.compile("\w .+\s*(\d+)\s*([\w\-]+)\s*(\d+)\s*\(\d+\)")
# Parses one HHresult body alignment line into one capturing group
BODY_ALIGNMENT_CONS = re.compile(
                                "\w .+\s*(\d+)\s*([\w\-~]+)\s*(\d+)\s*\(\d+\)")
# Parses one HHresult body alignment line into one capturing group
BODY_ALIGNMENT_MATCH = re.compile("\s+([\|\+\.\s*])")
# Parses one HHresult body alignment line into one capturing group
BODY_ALIGNMENT_CONF = re.compile("\w+\s+([\d\s]+)")
# ==============================================

HHMATCH_TABLE_ATTR = ["match_num", "hit_id", "probability", "expect",
                      "p_value", "score", "SS", "match_cols", "query_start",
                      "query_end", "hit_start", "hit_end", "hit_length"]
HHMATCH_BODY_ATTR = ["hit_id", "probability", "expect", "score", "match_cols",
                     "pid", "similarity", "sum_probs", "template_Neff"]

HHMATCH_ATTR = list(set(HHMATCH_TABLE_ATTR + HHMATCH_BODY_ATTR))


# ERROR HANDLING
# -----------------------------------------------------------------------------
HHRESULT_BASE_MESSAGE = ("Encountered improper formatting while parsing "
                         "HHResult file.\n")
HHRESULT_ERROR_MESSAGES = {"header": ("HHResult file header "
                                      "could not be recognized"),
                           "table": ("HHResult file data table "
                                     "could not be recognized"),
                           "body": ("HHResult file match result "
                                    "could not be recognized"),
                           "alignment": ("HHResult file match alignment "
                                         "was improperly formatted")}


class InitializationError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class HHResultFormatError(Exception):
    def __init__(self, e, line, line_num):
        line_traceback = "".join([f"at line {line_num}:\n",
                                  f"Line {line_num}> ", line])
        e = " ".join([e, line_traceback])
        super().__init__(e)


# HHRESULT DATA CLASSES
# -----------------------------------------------------------------------------
class HHResult:
    """Class for handling HHsuite result files.
    """
    def __init__(self, filepath):
        self.__filepath__ = filepath
        self.__initialized__ = False
        self.__lcounter = 1

        self.matches = None

        # Sets object attributes in preparatin for HHresult file headers
        for attr_name in HHRESULT_HEADERS.keys():
            setattr(self, attr_name, None)

    def parse_result(self):
        """Parses HHsuite result file.
        """
        with open(self.__filepath__, "r") as filehandle:
            self.__lcounter = 1
            self._parse_header(filehandle)
            match_index_map = self._parse_table(filehandle)
            self._parse_body(filehandle, match_index_map)

            self.matches = list(match_index_map.values())

    def _parse_header(self, filehandle):
        """Parses HHsuite result file header.

        :param filehandle: An open filehandle for a HHresult-formatted file
        """
        # For each HHresult file header line parse and store information
        for attr_name, reg_expr in HHRESULT_HEADERS.items():
            # Parses HHresult file header line and takes header_group
            header_line = self.__attempt_explicit_read(
                                        filehandle, reg_expr, ltype="header")
            header_split = re.split(reg_expr, header_line)

            # Sets attribute defined from HHRESULT_HEADERS with
            # expression capturing group
            setattr(self, attr_name, header_split[1])

        self.__lcounter += 1
        filehandle.readline()

    def _parse_table(self, filehandle):
        """Parses HHsuite result file table.

        :param filehandle: An open filehandle for a HHresult-formatted file
        :returns: Dictionary mapping HHresult match index to a HHMatch object
        :rtype: dict{}
        """
        self.__attempt_explicit_read(filehandle, MATCHES_HEADER, ltype="table")

        match_index_map = dict()
        # Reads matches until the tail end of a line that matches the
        # TABLE_DATA regular expression cannot be found
        while True:
            table_line, matches = self.__attempt_read_check(
                                            filehandle, TABLE_DATA)
            if not matches:
                break

            match_data_split = re.split(TABLE_DATA, table_line)
            match_id_split = re.split(TABLE_INDEX,
                                      match_data_split.pop(0))

            # Create match and use parsed information while ignoring
            # line end whitespace
            match = HHMatch(self.query_id)
            match.load_from_table_data(match_id_split[1:3] +
                                       match_data_split[0:11])

            match_index_map[match.match_num] = (match)

        return match_index_map

    def _parse_body(self, filehandle, match_index_map):
        self.__lcounter += 1
        body_line = filehandle.readline()
        while True:
            if re.match(BODY_INDEX, body_line) is None:
                break

            body_index_split = re.split(BODY_INDEX, body_line)

            body_line = self.__attempt_explicit_read(
                                    filehandle, BODY_NAME, ltype="body")
            body_name_split = re.split(BODY_NAME, body_line)

            body_line = self.__attempt_explicit_read(
                                    filehandle, BODY_DATA, ltype="body")
            body_data_split = re.split(BODY_DATA, body_line)

            match_index_map[body_index_split[1]].load_from_body_data(
                                        body_name_split[1:2] +
                                        body_data_split[1:-1])

            for _ in range(2):
                self.__lcounter += 1
                body_line = filehandle.readline()

            while True:
                if re.match(BODY_ALIGNMENT_SEQ, body_line) is not None:
                    alignment_data = self._parse_body_alignment(
                                            filehandle, body_line)
                    break

            for _ in range(2):
                self.__lcounter += 1
                body_line = filehandle.readline()

    def _parse_body_alignment(self, filehandle, body_line):
        query_seq_split = re.split(BODY_ALIGNMENT_SEQ, body_line)

        body_line = self.__attempt_explicit_read(
                                filehandle, BODY_ALIGNMENT_CONS,
                                ltype="alignment")
        query_cons_split = re.split(BODY_ALIGNMENT_CONS, body_line)

        body_line = self.__attempt_explicit_read(
                                filehandle, BODY_ALIGNMENT_MATCH,
                                ltype="alignment")
        match_lines_split = re.split(BODY_ALIGNMENT_MATCH, body_line)

        body_line = self.__attempt_explicit_read(
                                filehandle, BODY_ALIGNMENT_CONS,
                                ltype="alignment")
        target_cons_split = re.split(BODY_ALIGNMENT_CONS, body_line)

        body_line = self.__attempt_explicit_read(
                                filehandle, BODY_ALIGNMENT_SEQ,
                                ltype="alignment")
        target_seq_split = re.split(BODY_ALIGNMENT_SEQ, body_line)

        body_line = self.__attempt_explicit_read(
                                filehandle, BODY_ALIGNMENT_CONF,
                                ltype="alignment")
        match_confidence_split = re.split(BODY_ALIGNMENT_CONF, body_line)

    def __attempt_read_check(self, filehandle, regex):
        self.__lcounter += 1
        line = filehandle.readline()

        return (line, (re.search(regex, line) is not None))

    def __attempt_explicit_read(self, filehandle, regex, ltype=None):
        line = filehandle.readline()
        if re.search(regex, line) is None:
            e = HHRESULT_BASE_MESSAGE + HHRESULT_ERROR_MESSAGES.get(ltype, "")
            raise HHResultFormatError(e, line, self.__lcounter)

        self.__lcounter += 1

        return line

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

        self.seq_alignment = HHAlignment()
        self.cons_alignment = HHAlignment()

        # Sets object attributes in preparatin for HHresult file matches
        for i in range(len(HHMATCH_ATTR)):
            setattr(self, HHMATCH_ATTR[i], None)

    def load_from_table_data(self, table_data):
        # For each HHresult table match line column parse and store information
        for i in range(len(table_data)):
            setattr(self, HHMATCH_TABLE_ATTR[i], table_data[i])

    def load_from_body_data(self, body_data):
        # For each HHresult body match line column parse and store information
        for i in range(len(body_data)):
            setattr(self, HHMATCH_BODY_ATTR[i], body_data[i])


class HHAlignment:
    def __init__(self):
        self.alignment = None

    def extend_alignment(self, query_seq, target_seq):
        pass
