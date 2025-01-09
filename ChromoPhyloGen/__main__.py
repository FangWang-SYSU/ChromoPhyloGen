#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : __main__.py
@Author : XinWang
"""

import sys
from .cli_main import main
if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)
