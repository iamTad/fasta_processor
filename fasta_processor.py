"""
Class to access fasta data easily
"""
class Fasta():
    def __init__(self, fasta:str, fai:str):
        self.fasta = fasta # path to fasta file as we probe it
        self.fai = fai # path to fasta index file {contig: (startbyte, char/line, byte/line)}
        self.fai_dict = self.__parse_fai()

    def __parse_fai(self) -> dict[str, tuple[int, int, int, int]]:
        """
        Used in initialization. Loops through fai file and returns a dictionary of
        {contig: (startbyte, char per line, bytes per line, contig size)}
        """
        fai_dict = {}
        with open(self.fai, 'r') as fai_fs:
            for line in fai_fs:
                contig, contig_size, offset, base_per_line, byte_per_line = line.strip().split()
                contig_size_by_bytes = (int(contig_size) // int(base_per_line) * int(byte_per_line)) + (int(contig_size) % int(base_per_line))
                fai_dict[contig] = (int(offset), int(base_per_line), int(byte_per_line), int(contig_size))
        return fai_dict
   
    def get_sequence(self, contig:str, start:int, end:int):
        """
        Given a contig, sequence start (inclusive), and sequence end position (exclusive), 
        returns the corresponding sequence. This function allows for sequence extraction without
        storing the fasta file sequences.
        Start and end will be limited to 0 and the size of the contig
        """
        start = max(start, 0)
        end = min(end, self.fai_dict[contig][3])
        # ensure start less than end
        if start >= end:
            raise Exception("Start is not < end when extracting sequence from fasta")
        # estimate the number of lines we need to jump to after the header if multiline fasta
        lines = start // self.fai_dict[contig][1]
        lines_by_bytes = lines * self.fai_dict[contig][2]
        # get start byte
        start_leftover = start % self.fai_dict[contig][1] # how many bases left over
        start_byte = self.fai_dict[contig][0] + start_leftover + lines_by_bytes
        # get end byte
        seq_length = end - start
        seq_lines = seq_length // self.fai_dict[contig][1] # how many additional lines guaranteed to take up
        seq_lines_by_bytes = seq_lines * self.fai_dict[contig][2]
        seq_leftover = seq_length % self.fai_dict[contig][1] # how many extra bases left over in sequence after seq_lines
        seq_overflow = (start_leftover + seq_leftover) // self.fai_dict[contig][1] # will end leftover lead to new line
        # get calculate number of bytes to read after start byte
        bytes_to_read = seq_lines_by_bytes + seq_leftover + seq_overflow -1
        # extract sequence
        fasta = open(self.fasta, 'r')
        fasta.seek(start_byte) # move pointer to start byte
        sequence = fasta.read(bytes_to_read)
        sequence = ''.join(sequence.strip().split())
        return sequence
