class FastaFile:
    """class object for genbank flat
    file
    """

    def __init__(self, file):
        
        self.file = file
        self.body = list()

    def __enter__(self):
        """Database initialization 
        when initialized with 'With'
        """

        self.parse()
 

        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """'With' closurce
        """

        if exception_type is None:
            return True
        else:
            emessage = f"Error type: {exception_type}," \
                    f"Error Message: {exception_value}," \
                    f"Traceback: {traceback}"
            return emessage

    def parse(self):
        """Parse fasta file
        """

        try:
            header = ""
            sequence = ""
            header_count = 0

            # Loop file to read fasta file
            with open(self.file, "r") as f:
                for line in f:
                    if line.strip().startswith(">"):
                        header_count += 1
                        if header_count > 1:

                            # Save fasta format
                            self.body.append(
                                                {
                                                    "description":header,
                                                    "sequence":sequence.upper(),
                                                    "gc":(sequence.upper().count("G")+sequence.upper().count("C"))/len(sequence)
                                                }
                                            )

                            # Reset intermediate result
                            header_count = 1
                            header = ""
                            sequence = ""

                        header += line.strip().strip(">")
                    else:
                        sequence += line.strip("\n").strip("\t").strip()

                # Save fasta format
                if header:
                    self.body.append(
                                        {
                                            "description":header,
                                            "sequence":sequence.upper(),
                                            "gc":(sequence.upper().count("G")+sequence.upper().count("C"))/len(sequence)
                                        }
                                    )
        except Exception as e:
            msg = f"COULD NOT PARSE {self.file} | {type(e).__name__}, {e.args}, LINE AT: {e.traceback.tb_lineno}"
            print(msg)
    
    def collect(self):
        """Collect information
        """

        return self.body