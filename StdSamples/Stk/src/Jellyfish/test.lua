require('jam')
vm = jam.PythonVM()
jam.exec("print('hi')")
jam.execScript("test.py")
