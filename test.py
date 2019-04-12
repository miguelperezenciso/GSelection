#!/usr/bin/env python3
#-*- ecoding: utf-8 -*-

#Imports
from tkinter import *
import tkinter as tk
import tkinter.scrolledtext as tkst

mainWindow = tk.Tk()

#Config GUI
menubar = Menu(mainWindow)
mainWindow.config(menu=menubar)
ferramentas = Menu(menubar)
graph = Menu(menubar)

## Ferramentas
menubar.add_cascade(label='Ferramentas',font='12', menu=ferramentas)
menubar.add_separator()
ferramentas.add_command(label='Open Gen',font='18')
ferramentas.add_command(label='Resumo do Banco de Dados', font='18')
ferramentas.add_command(label='Gráficos',font='18')
ferramentas.add_command(label='Salvar',font='18')
ferramentas.add_separator()
ferramentas.add_command(label='Sair',font='18')

# Graph
menubar.add_cascade(label='Graphs',font='18', menu=graph)
graph.add_command(label='PCA',font='18')

#Output
saida = tkst.ScrolledText(master = mainWindow,wrap= WORD,width  = 20,height = 10)
saida.pack(padx=10, pady=10, fill=BOTH, expand=True)

if __name__== '__main__':
	mainWindow.title('GSelection')
	mainWindow.geometry('{}x{}'.format(mainWindow.winfo_screenwidth(),mainWindow.winfo_screenheight()))
	mainWindow.mainloop()