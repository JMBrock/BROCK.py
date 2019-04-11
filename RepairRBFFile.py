class brokenRBFile:

    def __init__(self,file):
        self.file = file

    def file_len(self):
        """
        This reads the length of a file
        """
        i = 0
        with open(self.file) as f:
            for i, l in enumerate(f):
                pass
        self.N_lines = i + 1
        print("-- Number of lines: %d" % self.N_lines )

    def repair_and_write(self,OutFile,stride):
        import ProgressBar
        with open('Rigid_Body.dat', 'r') as f:
            lines_ = f.readlines()
        lines = lines_[2:]

        f = open(OutFile, "wt" )
        self.N_Data = int(float(self.N_lines)/float(stride))

        for i in range(self.N_Data):
            for j in range(stride):
                values = lines[j+i*stride].strip('\n')
                f.write(values )
            f.write('\n')
            precent = float(i+1)/float(self.N_Data)*100.0
            ProgressBar.ProgressBar(i+1,self.N_Data,precent)

        print("")