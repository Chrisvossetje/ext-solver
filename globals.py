# Perform extra sanity checks on generated comodules
TEST = False

# Field to do algebra over
FIELD = 2

# Maximum filtration index our resolution will go to
FILTRATION_MAX = 20

# Maximum grade an element can occur in (wrt filtration index)
GRADE_LIMIT = (63, 0)

# This is necessary because you need some "space" around the grade_limit when creating an injection to a cofree comodule
ELEMENT_LIMIT = (GRADE_LIMIT[0] + 2, 0)