import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import linalg

node = hou.pwd()
geo = node.geometry()

# the node we are going to get our matrices from
cache_node = hou.node("../DEC_build_node")
cache_node2 = hou.node("../star0_inv")

# read heat
heat = np.array(geo.pointFloatAttribValues("heat"))
heat_prev = np.array(geo.pointFloatAttribValues("heat_prev"))
F = hou.intFrame()

# read time step
time_step = geo.attribValue("time_step")

# read matrices
LAP = cache_node2.cachedUserData("LAP")
star0 = cache_node.cachedUserData("star0")

# LHS = star0 - time_step*L
# RHS = star0.dot( heat )
# heat_new, tmp = linalg.cg(LHS,RHS)
for i in range(F):
    heat_prev = heat
    heat_new = time_step ** 2 * LAP.dot(heat) + star0.dot(2 * heat - heat_prev)

# update heat_prev attribute
# save to attribute
hou.Geometry.setPointFloatAttribValues(geo, "heat", heat_new)

# debug
# print "----"