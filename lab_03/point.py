class Point:

    def __init__(self, _x: float, _y: float) -> None:
        self.x = _x
        self.y = _y

    def __lt__(self, other):
        return self.x < other.x

    def __repr__(self):
        return "(x = {:.4f}; y = {:.4f})".format(self.x, self.y)


if __name__ == "__main__":
    print("This is package file")
