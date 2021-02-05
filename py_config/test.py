import Tkinter as tk
import tkMessageBox

class SampleApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.button = tk.Button(self, text="Push me", command=self.OnButton)
        self.button.pack()

    def OnButton(self):
        result = tkMessageBox.askokcancel(title="File already exists", 
                                       message="File already exists. Overwrite?")
        if result is True:
            print "User clicked Ok"
        else:
            print "User clicked Cancel"

if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()
