# The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO) 2023-2024.
from python_template_uv.init_utils.config import get_logger

logger = get_logger("example")


class ExampleClass:
    """
    This is a sample Python Template project class to demonstrate basic code and commenting style.
    """

    def __init__(self, place: str):
        """
        Creates a new ExampleClass to say hi from a ``place``
        :param place:
        """
        self.place = place

    def print_hi_from_place(self, name: str) -> str:
        """
        Prints hi to `name` from a `place`
        :param name: name to use
        return the printed string
        """
        str = f"Hi, {name}. From {self.place}"
        logger.info(str)
        return str

    @classmethod
    def print_hi(cls, name: str) -> None:
        """Says Hi to `name`"""
        logger.info(f"Hi, {name}")


def main() -> None:
    """Main entrypoint for the script"""
    # The class method can be called without instanciating the class
    ExampleClass.print_hi("You")

    # Instanciate the class at a location, then call the object's method
    ex = ExampleClass("Your Python Example Class")
    ex.print_hi_from_place("You")


if __name__ == "__main__":
    main()
