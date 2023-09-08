#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)

import sys

def loadingBar(count, total, task_part=None):
      """ Provides user with a loadingbar line. 
      See following:
            e.g. real)        041/400  [==                  ]             Subtask 793
            vars)     count/total [==                  ]             'task_part'
      Parameters
      ----------
      count : str, float or int
      Current task count. Easy to access throught 'enumerate()'
      total : str, float or int
            Maximal number of all tasks
      task_part : String | Optional
            If the task is divided in subtask and you want to keep track of
            your functions progress in detail pass your subtask in string format.
      Example
      -------
      array = np.linspace(1, 1000, 400)
      for p, i in enumerate(array):
      
      loadingBar(count=p, total=array.shape[0], task_part='Subtask')

      Returns
      -------
      stdout : Rewriteable String Output
            Generates a String Output for every of the progress steps
      """
      if task_part is None:
            task_part = ''
      percent = float(count + 1) / float(total) * 100
      size = 2

      sys.stdout.write("\r    "
                  + str(int(count + 1)).rjust(3, '0')
                  + "/" + str(int(total)).rjust(3, '0')
                  + ' [' + '=' * int(percent / 10) * size
                  + ' ' * (10 - int(percent / 10)) * size
                  + ']  %30s' % (task_part))
      if count + 1 == total:
            finish = '[done]'
            sys.stdout.write("\r    "
                        + str(int(count + 1)).rjust(3, '0')
                        + "/" + str(int(total)).rjust(3, '0')
                        + ' [' + '=' * int(percent / 10) * size
                        + ' ' * (10 - int(percent / 10)) * size
                        + ']  %30s\n' % (finish))

      return
