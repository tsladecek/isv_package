import pathlib


class Settings:
    def __init__(self):
        root_dir = str(pathlib.Path(__file__).parent.absolute())
        self.model_dir = root_dir + '/models/'
        self.data_dir = root_dir + '/data/'


settings = Settings()