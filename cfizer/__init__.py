from cfizer.utils import str_to_class, first_rest  # Allows use by all modules


class ConfigError(Exception):
    "Raised when something critical is missing from, or incorrectly defined in, config.yml. This should cause the program to exit immediately."
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class ConfigWarning(Exception):
    "Raised when something non-critical is missing from, or incorrectly defined in, config.yml. When this arises, it is intended that the program will be allowed to continue, but ending with a non-zero exit code."
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class VocabError(Exception):
    "Raised when something critical is missing from, or incorrectly defined in, vocabulary.yml. This should cause the program to exit immediately."
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
    

class VocabWarning(Exception):
    "Raised when something non-critical is missing from, or incorrectly defined in, vocabulary.yml. When this arises, it is intended that the program will be allowed to continue, but ending with a non-zero exit code."
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


