class DataClass():
    def _init__(self, data: dict = None):
        """
        Initialize the DataClass with a dictionary of data.
        
        :param data: Dictionary containing data to initialize the DataClass.
        """
        if data is None:
            data = {}
        self.data = data