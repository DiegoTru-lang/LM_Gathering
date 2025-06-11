class DataClass():
    def _init__(self, data: dict = None):
        """
        Initialize the DataClass with a dictionary of data.
        
        :param data: Dictionary containing data to initialize the DataClass.
        """
        if data is None:
            data = {}
        self.data = data

        """
        Debería leer todos los parámetros de un archivo Excel y almacenarlos. 
        Luego, en base a esos parámetros, y quizás algunos atributos, debería construir el resto de parámetros faltantes en función de esos parámetros
        También debería considerar transformar las unidades en las requeridas

        Finalmente, debería contener una clase Enum que me mapee "key" de mis parámetros a las variables en el modelo	
        """