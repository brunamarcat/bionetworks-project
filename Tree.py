
class Node:
    def _init_(self):
      self.left = None
      self.right = None
      self.cells = None
      self.type = None
      self.gens_split = None


    def to_json(self):
        # Create a dictionary to store the node's data
        data = {}

        # Check if attributes are not None before adding them to the dictionary
        if self.left is not None:
            data["left"] = self.left.to_json()
        if self.right is not None:
            data["right"] = self.right.to_json()
        if self.cells is not None:
            data["cells"] = self.cells
        if self.type is not None:
            data["type"] = self.type
        if self.gens_split is not None:
            data["gens_split"] = self.gens_split
        return data

    def from_json(self, data):
        # Check if keys exist in the data dictionary before accessing them
        self.left = Node() if "left" in data else None
        if self.left:
            self.left.from_json(data["left"])
        self.right = Node() if "right" in data else None
        if self.right:
            self.right.from_json(data["right"])
        self.cells = data.get("cells")
        self.type = data.get("type")
        self.gens_split = data.get("gens_split")

