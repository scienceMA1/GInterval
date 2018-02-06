class Parser(object):
    """
    `Parser` is the basic super class for the creation of
    GInterval instance from file. Any sub-class should
    extend from this class and constructed by manually or
    :class:`ParserFactory`.

    Args:
        path (str): The path of the file.

    """
    def __init__(self, path):
        self.path = path
