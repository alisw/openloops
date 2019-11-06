#!/usr/bin/env python
# -*- coding: utf-8 -*-

#!******************************************************************************!
#! Copyright (C) 2014-2019 OpenLoops Collaboration. For authors see authors.txt !
#!                                                                              !
#! This file is part of OpenLoops.                                              !
#!                                                                              !
#! OpenLoops is free software: you can redistribute it and/or modify            !
#! it under the terms of the GNU General Public License as published by         !
#! the Free Software Foundation, either version 3 of the License, or            !
#! (at your option) any later version.                                          !
#!                                                                              !
#! OpenLoops is distributed in the hope that it will be useful,                 !
#! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
#! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                !
#! GNU General Public License for more details.                                 !
#!                                                                              !
#! You should have received a copy of the GNU General Public License            !
#! along with OpenLoops.  If not, see <http://www.gnu.org/licenses/>.           !
#!******************************************************************************!

import multiprocessing


class WorkerProcess(multiprocessing.Process):
    """A worker process to evaluate matrix elements."""

    def __init__(self, processes, options, swap_options, in_queue, out_queue):
        """Create a worker subprocess to evaluate matrix elements in a
        subprocess.
        processes -- [(proc, amptype), ...] processes as strings with amptypes
        options -- [(key,val), ...] to initialise with set_parameter()
        swap_options -- [(key, [val1, val2, ...]), ...] options for
                        multiple evaluations with different parameter values
        in_queue -- multiprocessing.Queue() to receive
                    (n, (process, amptype), PhaseSpacePoint) objects
        out_queue -- multiprocessing.Queue() to put MatrixElement objects
        """
        super(WorkerProcess, self).__init__()
        self.processes = processes
        self.options = options
        self.swap_options = swap_options
        self.in_queue = in_queue
        self.out_queue = out_queue

    def run(self):
        """Start the worker subprocess, set parameters and register processes.
        Run loop: receive (data_id, (process_str, amptype), phase space point)
        from in_queue and evaluate the matrix element (multiply if swap
        optionsare set). Push the result(s) as a list of matrix elements
        to out_queue."""
        proc_handles = {}
        import openloops
        for key, val in self.options:
            openloops.set_parameter(key, val)
        for procamp in self.processes:
            proc_handles[procamp] = openloops.Process(*procamp)
        while True:
            data_id, procamp, psp = self.in_queue.get()
            mes = []
            if self.swap_options:
                # multiple evaluation with different parameters
                for sopts in self.swap_options:
                    for key, val in sopts:
                        openloops.set_parameter(key, val)
                    mes.append(proc_handles[procamp].evaluate(psp))
            else:
                # single evaluation
                mes.append(proc_handles[procamp].evaluate(psp))
            self.out_queue.put((data_id, mes))


class Pool(object):
    """A pool of worker processes to evaluate matrix elements in parallel."""

    def __init__(self, processes, workers=0, options=[], swap_options=[]):
        """Create worker processes in the pool and start them as daemons."""
        if not workers:
            self.n_workers = multiprocessing.cpu_count()
        else:
            self.n_workers = workers
        self.in_queue = multiprocessing.Queue()
        self.out_queue = multiprocessing.Queue()
        # number of submitted points
        self.n_put = 0
        # number of received points
        self.n_get = 0
        self.res_buffer = {}
        # start worker processes
        for n in range(self.n_workers):
            worker = WorkerProcess(processes, options, swap_options,
                                   self.in_queue, self.out_queue)
            # die with the main process
            worker.daemon = True
            worker.start()

    def put(self, procamp, psp):
        """Submit a ((process, amptype), phase_space_point)."""
        self.in_queue.put((self.n_put, procamp, psp))
        self.n_put += 1

    def get(self):
        """Receive a MatrixElement.
        Matrix elements are received in the order in which the corresponding
        phase space points were submitted."""
        me = self.res_buffer.pop(self.n_get, None)
        if me:
            self.n_get += 1
            return me
        while True:
            data_id, me = self.out_queue.get()
            if data_id == self.n_get:
                self.n_get += 1
                return me
            else:
                self.res_buffer[data_id] = me
