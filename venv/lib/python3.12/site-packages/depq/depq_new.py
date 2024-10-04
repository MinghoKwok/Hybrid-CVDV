from collections import deque


class DoublyLinkedListNode:
    def __init__(self, key, data):
        self.key = key
        self.data = data
        self.prev = None
        self.next = None

    def __repr__(self):
        return 'Node({}, {})'.format(self.key, self.data)

    def __str__(self):
        return self.__repr__()


class DoublyLinkedList:
    def __init__(self):
        self.first = None
        self.last = None
        self.size = 0

    def append(self, key, data):

        if not self.last:
            self.first = self.last = DoublyLinkedListNode(key, data)
        else:
            prev_last = self.last
            last = self.last = prev_last.next = DoublyLinkedListNode(key, data)
            last.prev = prev_last
            last.next = self.first

        self.size += 1

    def appendleft(self, key, data):

        if not self.first:
            self.first = self.last = DoublyLinkedListNode(key, data)
        else:
            prev_first = self.first
            first = self.first = prev_first.prev = DoublyLinkedListNode(key, data)
            first.prev = self.last
            first.next = prev_first

        self.size += 1

    def __iter__(self):
        node = self.first

        for _ in range(self.size):
            yield node
            node = node.next

        raise StopIteration

    def __repr__(self):
        return ', '.join(str(s) for s in iter(self))

    def __str__(self):
        return self.__repr__()


class DEPQ:
    def __init__(self):
        pass

    def popfirst(self):
        pass





