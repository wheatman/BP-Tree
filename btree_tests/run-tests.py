import os
import tarfile
from subprocess import call
import sys

for i in range(1, 9):
  n = 10**i
  call(["./basic", str(n), str(5)])
