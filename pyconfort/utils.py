"""
This module contains some classes and functions that are used from other modules
"""

#class for logging
class Logger:
    """
    Class that wraps a file object to abstract the logging.
    """
    # Class Logger to writargs.input.split('.')[0] output to a file
    def __init__(self, filein, append, suffix='dat'):
        self.log = open(f'{filein}_{append}.{suffix}', 'w')

    def write(self, message):
        """
        Appends a newline character to the message and writes it into the file.
   
        Parameters
        ----------
        message : str
           text to be written in the log file.
        """
        self.log.write(f'{message}\n')

    def fatal(self, message):
        """
        Writes the message to the file. Closes the file and raises an error exit
   
        Parameters
        ----------
        message : str
           text to be written in the log file.
        """
        self.write(message)
        self.finalize()
        raise SystemExit(1)

    def finalize(self):
        """ Closes the file """
        self.log.close()