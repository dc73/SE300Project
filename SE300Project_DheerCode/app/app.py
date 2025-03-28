class App:
    def __init__(self):
        from gui import MainApp 
        self.gui = MainApp()

    def run(self):
        self.gui.mainloop()


if __name__ == "__main__":
    app = App()
    app.run()
