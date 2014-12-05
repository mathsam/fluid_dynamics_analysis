def cmp_list(list1, list2):
    """Lists are compared according to the first element first. If equal, then
    use the second element, and so on. For example
    [1, 3] > [0, 0]
    [1, 3] > [1, 2]
    [1.]   > [0, 2]
    [1,2]  > [1,]
    non-empty list > []

    return -1 if list1 <  list2
    return 0  if list1 == list2
    return 1  if list1 >  list2

    Note: expect items in list to be strings representing integers, so convert
    strings into integers and then compare.
    """
    if list1 == list2:
        return 0
    for item1, item2 in zip(list1, list2):
        if(int(item1) > int(item2)): return 1
        if(int(item1) < int(item2)): return -1
    return 1 if len(list1)>len(list2) else -1
 
class ComparableList(object):
    """Turn list into comparable object."""
    def __init__(self, obj, *args):
	self.obj = obj
    def __lt__(self, other):
	return cmp_list(self.obj, other.obj) < 0
    def __gt__(self, other):
	return cmp_list(self.obj, other.obj) > 0
    def __eq__(self, other):
	return cmp_list(self.obj, other.obj) == 0
    def __le__(self, other):
	return cmp_list(self.obj, other.obj) <= 0
    def __ge__(self, other):
	return cmp_list(self.obj, other.obj) >= 0
    def __ne__(self, other):
	return cmp_list(self.obj, other.obj) != 0
