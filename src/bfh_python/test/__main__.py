"""

There are 29 modules in the main level of bfh_python.  The following
have no tests: [__init__, autocompleteda, experimental, data, tests]

"""

import importlib
import unittest
import sys
import time

# The modules are ordered so that each module is tested only after all
# of the modules it imports.

module_names = ['utility',
                'latex',
                'algebra',
                'linalg',
                'grading',
                'pmc',
                'signs',
                'minusalg',
                'localpmc',
                'hdiagram',
                'identityaa',
                'dstructure',
                'ddstructure',
                'dastructure',
                'extendbyid',
                'digraph',
                'arcslide',
                'cobordism',
                'dehntwist',
                'arcslideda',
                'dehntwistda',
                'cobordismda',
                'involutive',
                'braid']

if __name__ == '__main__':
    start_time = time.time()
    modules_done, tests_run, errors, failures = 0, 0, 0, 0
    for name in module_names:
        print(f'*** Starting tests of module {name} ***\n')
        module = importlib.import_module('.' + name + 'test', package=__package__)
        result = unittest.main(module, exit=False).result
        modules_done += 1
        tests_run += result.testsRun
        errors += len(result.errors)
        failures += len(result.failures)
        time_elapsed = time.time() - start_time
        print(f'\nCumulative: modules={modules_done} run={tests_run} errors={errors} failures={failures} time={time_elapsed:0.3f}s\n')

    sys.exit(errors + failures)

    
