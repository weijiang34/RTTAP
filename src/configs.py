import os

class QCDir:
    def __init__(self, output_dir):
        self._path = os.path.join(output_dir, 'quality_control')
        self.clean_fq = os.path.join(self._path, '{fileHeader}.clean.fq.gz')
    
    @property
    def path(self):
        return self._path