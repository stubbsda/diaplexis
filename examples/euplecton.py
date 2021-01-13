#!/usr/bin/env python3

# Use the command 
# sudo pip3 install tkinter 
# to install the Tkinter module for Python 3.x.
import tkinter
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
import os
import math
import lxml

# One point to note is that this GUI must be run on a monitor with
# resolution at least 800 x 675 for proper viewing of the complete interface. 

class euplecton:
    def __init__(self,master=None):
        self.master = master
        self.master.title('Euplecton')
        self.master.geometry('730x850')
        self.master.resizable(0,0)

        self.sheet_dynamics = tkinter.BooleanVar()
        self.superposable = tkinter.BooleanVar()
        self.compressible = tkinter.BooleanVar()
        self.permutable = tkinter.BooleanVar()
        self.perturbg = tkinter.BooleanVar()
        self.perturbt = tkinter.BooleanVar()
        self.perturbe = tkinter.BooleanVar()
        self.dim_uniformity = tkinter.BooleanVar()
        self.relational = tkinter.BooleanVar()
        self.high_memory_footprint = tkinter.BooleanVar()
        self.cg_refinement = tkinter.BooleanVar()

        self.initial_state = tkinter.StringVar()
        self.hyphansis = tkinter.StringVar()
        self.homology_method = tkinter.StringVar()
        self.homology_field = tkinter.StringVar()
        self.signature = tkinter.StringVar()
        self.solver_type = tkinter.StringVar()
        self.parameter_filename = tkinter.StringVar()
        self.input_filename = tkinter.StringVar()
        self.hyphansis_score = tkinter.StringVar()
        self.int_engine = tkinter.StringVar()

        self.initial_events = tkinter.IntVar()
        self.max_iterations = tkinter.IntVar()
        self.initial_dimension = tkinter.IntVar()
        self.initial_sheets = tkinter.IntVar()
        self.chkpt_frequency = tkinter.IntVar()
        self.background_dim = tkinter.IntVar()
        self.random_seed = tkinter.IntVar()
        self.solver_iterations = tkinter.IntVar()
        self.max_generations = tkinter.IntVar()
        self.pool_size = tkinter.IntVar()
        self.max_jousts = tkinter.IntVar()
        self.thermal_sweeps = tkinter.IntVar()
        self.annealing_steps = tkinter.IntVar()
        self.max_int_steps = tkinter.IntVar()
        self.max_cg_steps = tkinter.IntVar()
        self.max_ls_steps = tkinter.IntVar()

        self.edge_probability = tkinter.DoubleVar()
        self.abnormality_threshold = tkinter.DoubleVar()
        self.superposition_threshold = tkinter.DoubleVar()
        self.parity_probability = tkinter.DoubleVar()
        self.geometry_threshold = tkinter.DoubleVar()
        self.thermal_variance = tkinter.DoubleVar()
        self.thermalization_criterion = tkinter.DoubleVar()
        self.step_size = tkinter.DoubleVar()
        self.damping = tkinter.DoubleVar()
        self.spring = tkinter.DoubleVar()
        self.repulsion = tkinter.DoubleVar()
        self.edge_flexibility = tkinter.DoubleVar()
        self.reflection = tkinter.DoubleVar()
        self.contraction = tkinter.DoubleVar()
        self.expansion = tkinter.DoubleVar()
        self.shrinkage = tkinter.DoubleVar()

        self.clear_parameters()

        global_group = tkinter.LabelFrame(self.master,text="Global",padx=5,pady=5)
        geometry_group = tkinter.LabelFrame(self.master,text="Geometry",padx=5,pady=5)
        button_group = tkinter.LabelFrame(self.master,text="",padx=5,pady=5,bd=0)

        global_group.grid(row=0,column=0,padx=5,pady=5,sticky="")
        geometry_group.grid(row=1,column=0,sticky="")
        button_group.grid(row=2,column=0,sticky="")

        initial_states = ['Cartesian','Singleton','Monoplex','Random']
        hyphansis_types = ['Dynamic','Musical']
        homology_methods = ['GAP','Native']
        homology_fields = ['GF2','NTL::ZZ','INT']
        signatures = ['Euclidean','Lorentzian']
        solvers = ['Minimal','Evolutionary','Annealing','Mechanical','Simplex']
        engines = ['Euler','RK4']

        self.label1 = tkinter.Label(global_group,text='Number of Initial Events:',wraplength=250,justify=tkinter.LEFT)
        label2 = tkinter.Label(global_group,text='Maximum Number of Relaxation Steps:',wraplength=250,justify=tkinter.LEFT)
        self.label3 = tkinter.Label(global_group,text='Initial Dimension:',wraplength=250,justify=tkinter.LEFT)
        label4 = tkinter.Label(global_group,text='Initial State:',wraplength=250,justify=tkinter.LEFT)
        initial_state = tkinter.OptionMenu(global_group,self.initial_state,*initial_states,command=self.istate_change)
        self.sheet_check = tkinter.Checkbutton(global_group,text='Sheet Dynamics',variable=self.sheet_dynamics)
        self.label6 = tkinter.Label(global_group,text='Number of Initial Sheets:',wraplength=250,justify=tkinter.LEFT)
        superposition_check = tkinter.Checkbutton(global_group,text='Superposable',variable=self.superposable)
        compression_check = tkinter.Checkbutton(global_group,text='Compressible',variable=self.compressible)
        permutation_check = tkinter.Checkbutton(global_group,text='Permutable',variable=self.permutable)
        label9 = tkinter.Label(global_group,text='Checkpoint Frequency:',wraplength=250,justify=tkinter.LEFT)
        perturb_geometry = tkinter.Checkbutton(global_group,text='Perturb Geometry',variable=self.perturbg)
        perturb_topology = tkinter.Checkbutton(global_group,text='Perturb Topology',variable=self.perturbt)
        perturb_energy = tkinter.Checkbutton(global_group,text='Perturb Energy',variable=self.perturbe)
        label11 = tkinter.Label(global_group,text='Hyphansis Type:',wraplength=250,justify=tkinter.LEFT)
        hyphansis_type = tkinter.OptionMenu(global_group,self.hyphansis,*hyphansis_types,command=self.htype_change)
        self.label12 = tkinter.Label(global_group,text='Edge Probability:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        label13 = tkinter.Label(global_group,text='Abnormality Threshold:',wraplength=250,justify=tkinter.LEFT)
        self.label14 = tkinter.Label(global_group,text='Parity Mutation Probability:',wraplength=250,justify=tkinter.LEFT)
        label15 = tkinter.Label(global_group,text='Superposition Threshold:',wraplength=250,justify=tkinter.LEFT)
        label17 = tkinter.Label(global_group,text='Homology Method:',wraplength=250,justify=tkinter.LEFT)
        homology_method = tkinter.OptionMenu(global_group,self.homology_method,*homology_methods)
        label18 = tkinter.Label(global_group,text='Homology Field:',wraplength=250,justify=tkinter.LEFT)
        homology_field = tkinter.OptionMenu(global_group,self.homology_field,*homology_fields)
        label20 = tkinter.Label(global_group,text='Random Number Generator Seed:',wraplength=250,justify=tkinter.LEFT)
        label21 = tkinter.Label(global_group,text='Spacetime Signature:',wraplength=250,justify=tkinter.LEFT)
        signature = tkinter.OptionMenu(global_group,self.signature,*signatures)
        footprint_check = tkinter.Checkbutton(global_group,text='Large Memory Footprint',variable=self.high_memory_footprint)
        uniformity_check = tkinter.Checkbutton(global_group,text='Dimensional Uniformity',variable=self.dim_uniformity)
        relational_check = tkinter.Checkbutton(global_group,text='Relational Geometry',variable=self.relational)
        label25 = tkinter.Label(global_group,text='Background Dimension:',wraplength=250,justify=tkinter.LEFT)
        self.label26 = tkinter.Label(global_group,text='Hyphansis Score:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)

        label30 = tkinter.Label(geometry_group,text='Solver Type:',wraplength=250,justify=tkinter.LEFT)
        solver_type = tkinter.OptionMenu(geometry_group,self.solver_type,*solvers,command=self.gsolver_change)
        label31 = tkinter.Label(geometry_group,text='Geometry Tolerance:',wraplength=250,justify=tkinter.LEFT)
        self.label32 = tkinter.Label(geometry_group,text='Solver Iterations:',wraplength=250,justify=tkinter.LEFT)
        self.label33 = tkinter.Label(geometry_group,text='Maximum Number of Generations:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label34 = tkinter.Label(geometry_group,text='Population Pool Size:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label35 = tkinter.Label(geometry_group,text='Maximum Number of Jousts:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label36 = tkinter.Label(geometry_group,text='Thermal Variance:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label37 = tkinter.Label(geometry_group,text='Number of Thermal Sweeps:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label38 = tkinter.Label(geometry_group,text='Thermalization Criterion:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label39 = tkinter.Label(geometry_group,text='Number of Annealing Steps:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label40 = tkinter.Label(geometry_group,text='Integration Engine:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.engine = tkinter.OptionMenu(geometry_group,self.int_engine,*engines)
        self.engine.config(state=tkinter.DISABLED)
        self.label41 = tkinter.Label(geometry_group,text='Step Size:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label42 = tkinter.Label(geometry_group,text='Maximum Number of Integration Steps:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label43 = tkinter.Label(geometry_group,text='Damping Constant:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label44 = tkinter.Label(geometry_group,text='Spring Constant:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label45 = tkinter.Label(geometry_group,text='Repulsion Constant:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.cg_check = tkinter.Checkbutton(geometry_group,text='Conjugate Gradient Refinement',variable=self.cg_refinement,state=tkinter.DISABLED)
        self.label46 = tkinter.Label(geometry_group,text='Maximum Number of Conjugate Gradient Steps:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label47 = tkinter.Label(geometry_group,text='Maximum Number of Line Solver Steps:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label48 = tkinter.Label(geometry_group,text='Edge Flexibility Threshold:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label49 = tkinter.Label(geometry_group,text='Simplex Reflection Coefficient:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label50 = tkinter.Label(geometry_group,text='Simplex Expansion Coefficient:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label51 = tkinter.Label(geometry_group,text='Simplex Contraction Coefficient:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
        self.label52 = tkinter.Label(geometry_group,text='Simplex Shrinkage Coefficient:',wraplength=250,justify=tkinter.LEFT,state=tkinter.DISABLED)
       
        self.entry1 = tkinter.Entry(global_group,width=7,textvariable=self.initial_events)
        entry2 = tkinter.Entry(global_group,width=7,textvariable=self.max_iterations)
        self.entry3 = tkinter.Entry(global_group,width=7,textvariable=self.initial_dimension)
        self.entry6 = tkinter.Entry(global_group,width=7,textvariable=self.initial_sheets)
        entry9 = tkinter.Entry(global_group,width=7,textvariable=self.chkpt_frequency)
        self.entry12 = tkinter.Entry(global_group,width=7,textvariable=self.edge_probability,state=tkinter.DISABLED)
        entry13 = tkinter.Entry(global_group,width=7,textvariable=self.abnormality_threshold)
        self.entry14 = tkinter.Entry(global_group,width=7,textvariable=self.parity_probability)
        entry15 = tkinter.Entry(global_group,width=7,textvariable=self.superposition_threshold)
        entry20 = tkinter.Entry(global_group,width=7,textvariable=self.random_seed)
        entry25 = tkinter.Entry(global_group,width=7,textvariable=self.background_dim)
        self.entry26 = tkinter.Entry(global_group,width=18,textvariable=self.hyphansis_score,state=tkinter.DISABLED)
        entry31 = tkinter.Entry(geometry_group,width=7,textvariable=self.geometry_threshold)
        self.entry32 = tkinter.Entry(geometry_group,width=7,textvariable=self.solver_iterations)
        self.entry33 = tkinter.Entry(geometry_group,width=7,textvariable=self.max_generations,state=tkinter.DISABLED)
        self.entry34 = tkinter.Entry(geometry_group,width=7,textvariable=self.pool_size,state=tkinter.DISABLED)
        self.entry35 = tkinter.Entry(geometry_group,width=7,textvariable=self.max_jousts,state=tkinter.DISABLED)
        self.entry36 = tkinter.Entry(geometry_group,width=7,textvariable=self.thermal_variance,state=tkinter.DISABLED)
        self.entry37 = tkinter.Entry(geometry_group,width=7,textvariable=self.thermal_sweeps,state=tkinter.DISABLED)
        self.entry38 = tkinter.Entry(geometry_group,width=7,textvariable=self.thermalization_criterion,state=tkinter.DISABLED)
        self.entry39 = tkinter.Entry(geometry_group,width=7,textvariable=self.annealing_steps,state=tkinter.DISABLED)
        self.entry41 = tkinter.Entry(geometry_group,width=7,textvariable=self.step_size,state=tkinter.DISABLED)
        self.entry42 = tkinter.Entry(geometry_group,width=7,textvariable=self.max_int_steps,state=tkinter.DISABLED)
        self.entry43 = tkinter.Entry(geometry_group,width=7,textvariable=self.damping,state=tkinter.DISABLED)
        self.entry44 = tkinter.Entry(geometry_group,width=7,textvariable=self.spring,state=tkinter.DISABLED)
        self.entry45 = tkinter.Entry(geometry_group,width=7,textvariable=self.repulsion,state=tkinter.DISABLED)
        self.entry46 = tkinter.Entry(geometry_group,width=7,textvariable=self.max_cg_steps,state=tkinter.DISABLED)
        self.entry47 = tkinter.Entry(geometry_group,width=7,textvariable=self.max_ls_steps,state=tkinter.DISABLED)
        self.entry48 = tkinter.Entry(geometry_group,width=7,textvariable=self.edge_flexibility,state=tkinter.DISABLED)
        self.entry49 = tkinter.Entry(geometry_group,width=7,textvariable=self.reflection,state=tkinter.DISABLED)
        self.entry50 = tkinter.Entry(geometry_group,width=7,textvariable=self.expansion,state=tkinter.DISABLED)
        self.entry51 = tkinter.Entry(geometry_group,width=7,textvariable=self.contraction,state=tkinter.DISABLED)
        self.entry52 = tkinter.Entry(geometry_group,width=7,textvariable=self.shrinkage,state=tkinter.DISABLED)

        label19 = tkinter.Label(button_group,text='Parameter Filename:',wraplength=250,justify=tkinter.LEFT)
        entry19 = tkinter.Entry(button_group,width=18,textvariable=self.parameter_filename)
        button1 = tkinter.Button(button_group,text='Save Parameters',command=self.save_parameters)
        button2 = tkinter.Button(button_group,text='Exit',command=root.quit)
        button3 = tkinter.Button(button_group,text='Load Parameters',command=self.load_parameters)
        button4 = tkinter.Button(button_group,text='Restore Default Values',command=self.clear_parameters)

        self.label1.grid(row=0,column=0,sticky=tkinter.W)
        self.entry1.grid(row=0,column=1,sticky=tkinter.W)
        label2.grid(row=1,column=0,sticky=tkinter.W)
        entry2.grid(row=1,column=1,sticky=tkinter.W)
        self.label3.grid(row=2,column=0,sticky=tkinter.W)
        self.entry3.grid(row=2,column=1,sticky=tkinter.W)
        self.label6.grid(row=3,column=0,sticky=tkinter.W)
        self.entry6.grid(row=3,column=1,sticky=tkinter.W)
        self.sheet_check.grid(row=4,column=0,sticky=tkinter.W)
        label4.grid(row=5,column=0,sticky=tkinter.W)
        initial_state.grid(row=5,column=1,sticky=tkinter.W)
        label13.grid(row=6,column=0,sticky=tkinter.W)
        entry13.grid(row=6,column=1,sticky=tkinter.W)
        label17.grid(row=7,column=0,sticky=tkinter.W)
        homology_method.grid(row=7,column=1,sticky=tkinter.W)
        label18.grid(row=8,column=0,sticky=tkinter.W)
        homology_field.grid(row=8,column=1,sticky=tkinter.W)
        label20.grid(row=9,column=0,sticky=tkinter.W)
        entry20.grid(row=9,column=1,sticky=tkinter.W)
        label9.grid(row=10,column=0,sticky=tkinter.W)
        entry9.grid(row=10,column=1,sticky=tkinter.W)
        self.label26.grid(row=11,column=0,sticky=tkinter.W)
        self.entry26.grid(row=11,column=1,sticky=tkinter.W)
        superposition_check.grid(row=0,column=3,sticky=tkinter.W)
        compression_check.grid(row=1,column=3,sticky=tkinter.W)
        permutation_check.grid(row=2,column=3,sticky=tkinter.W)
        label11.grid(row=3,column=3,sticky=tkinter.W)
        hyphansis_type.grid(row=3,column=4,sticky=tkinter.W)
        self.label12.grid(row=4,column=3,sticky=tkinter.W) 
        self.entry12.grid(row=4,column=4,sticky=tkinter.W)
        self.label14.grid(row=5,column=3,sticky=tkinter.W)
        self.entry14.grid(row=5,column=4,sticky=tkinter.W)
        label15.grid(row=6,column=3,sticky=tkinter.W)
        entry15.grid(row=6,column=4,sticky=tkinter.W)
        label21.grid(row=7,column=3,sticky=tkinter.W)
        signature.grid(row=7,column=4,sticky=tkinter.W)
        relational_check.grid(row=8,column=3,sticky=tkinter.W)
        uniformity_check.grid(row=9,column=3,sticky=tkinter.W)
        footprint_check.grid(row=10,column=3,sticky=tkinter.W)
        label25.grid(row=11,column=3,sticky=tkinter.W)
        entry25.grid(row=11,column=4,sticky=tkinter.W)

        label30.grid(row=0,column=0,sticky=tkinter.W)
        solver_type.grid(row=0,column=1,sticky=tkinter.W)
        label31.grid(row=1,column=0,sticky=tkinter.W)
        entry31.grid(row=1,column=1,sticky=tkinter.W)
        self.label40.grid(row=2,column=0,sticky=tkinter.W)
        self.engine.grid(row=2,column=1,sticky=tkinter.W)
        self.label41.grid(row=3,column=0,sticky=tkinter.W)
        self.entry41.grid(row=3,column=1,sticky=tkinter.W)
        self.label42.grid(row=4,column=0,sticky=tkinter.W)
        self.entry42.grid(row=4,column=1,sticky=tkinter.W)
        self.label43.grid(row=5,column=0,sticky=tkinter.W)
        self.entry43.grid(row=5,column=1,sticky=tkinter.W)
        self.label44.grid(row=6,column=0,sticky=tkinter.W)
        self.entry44.grid(row=6,column=1,sticky=tkinter.W)
        self.label45.grid(row=7,column=0,sticky=tkinter.W)
        self.entry45.grid(row=7,column=1,sticky=tkinter.W)
        self.cg_check.grid(row=8,column=0,sticky=tkinter.W)
        self.label46.grid(row=9,column=0,sticky=tkinter.W)
        self.entry46.grid(row=9,column=1,sticky=tkinter.W)
        self.label47.grid(row=10,column=0,sticky=tkinter.W)
        self.entry47.grid(row=10,column=1,sticky=tkinter.W)
        self.label48.grid(row=11,column=0,sticky=tkinter.W)
        self.entry48.grid(row=11,column=1,sticky=tkinter.W)

        self.label32.grid(row=0,column=3,sticky=tkinter.W)
        self.entry32.grid(row=0,column=4,sticky=tkinter.W)
        self.label33.grid(row=1,column=3,sticky=tkinter.W)
        self.entry33.grid(row=1,column=4,sticky=tkinter.W)
        self.label34.grid(row=2,column=3,sticky=tkinter.W)
        self.entry34.grid(row=2,column=4,sticky=tkinter.W)
        self.label35.grid(row=3,column=3,sticky=tkinter.W)
        self.entry35.grid(row=3,column=4,sticky=tkinter.W)
        self.label36.grid(row=4,column=3,sticky=tkinter.W)
        self.entry36.grid(row=4,column=4,sticky=tkinter.W)
        self.label37.grid(row=5,column=3,sticky=tkinter.W)
        self.entry37.grid(row=5,column=4,sticky=tkinter.W)
        self.label38.grid(row=6,column=3,sticky=tkinter.W)
        self.entry38.grid(row=6,column=4,sticky=tkinter.W)
        self.label39.grid(row=7,column=3,sticky=tkinter.W)
        self.entry39.grid(row=7,column=4,sticky=tkinter.W)
        self.label49.grid(row=8,column=3,sticky=tkinter.W)
        self.entry49.grid(row=8,column=4,sticky=tkinter.W)
        self.label50.grid(row=9,column=3,sticky=tkinter.W)
        self.entry50.grid(row=9,column=4,sticky=tkinter.W)
        self.label51.grid(row=10,column=3,sticky=tkinter.W)
        self.entry51.grid(row=10,column=4,sticky=tkinter.W)
        self.label52.grid(row=11,column=3,sticky=tkinter.W)
        self.entry52.grid(row=11,column=4,sticky=tkinter.W)

        button1.grid(row=0,column=0)
        button4.grid(row=0,column=1)
        button3.grid(row=1,column=0)
        button2.grid(row=1,column=1)
        label19.grid(row=2,column=0)
        entry19.grid(row=2,column=1)

    def htype_change(self,*args):
        if self.hyphansis.get() == 'Dynamic':
           self.sheet_check.config(state=tkinter.NORMAL)
           self.label6.config(state=tkinter.NORMAL)
           self.entry6.config(state=tkinter.NORMAL)
           self.label14.config(state=tkinter.NORMAL)
           self.entry14.config(state=tkinter.NORMAL)
           self.label26.config(state=tkinter.DISABLED)
           self.entry26.config(state=tkinter.DISABLED)
        else:
           self.sheet_check.config(state=tkinter.DISABLED)
           self.label6.config(state=tkinter.DISABLED)
           self.entry6.config(state=tkinter.DISABLED)
           self.label14.config(state=tkinter.DISABLED)
           self.entry14.config(state=tkinter.DISABLED)
           self.label26.config(state=tkinter.NORMAL)
           self.entry26.config(state=tkinter.NORMAL)

    def istate_change(self,*args):
        if self.initial_state.get() == 'Cartesian':
           self.label1.config(state=tkinter.NORMAL)
           self.entry1.config(state=tkinter.NORMAL)
           self.label3.config(state=tkinter.NORMAL)
           self.entry3.config(state=tkinter.NORMAL)
           self.label12.config(state=tkinter.DISABLED)
           self.entry12.config(state=tkinter.DISABLED)
        elif self.initial_state.get() == 'Singleton':
           self.label1.config(state=tkinter.DISABLED)
           self.entry1.config(state=tkinter.DISABLED)
           self.label3.config(state=tkinter.DISABLED)
           self.entry3.config(state=tkinter.DISABLED)
           self.label12.config(state=tkinter.DISABLED)
           self.entry12.config(state=tkinter.DISABLED)
        elif self.initial_state.get() == 'Monoplex':
           self.label1.config(state=tkinter.DISABLED)
           self.entry1.config(state=tkinter.DISABLED)
           self.label3.config(state=tkinter.NORMAL)
           self.entry3.config(state=tkinter.NORMAL)
           self.label12.config(state=tkinter.DISABLED)
           self.entry12.config(state=tkinter.DISABLED)
        elif self.initial_state.get() == 'Random':
           self.label1.config(state=tkinter.NORMAL)
           self.entry1.config(state=tkinter.NORMAL)
           self.label3.config(state=tkinter.NORMAL)
           self.entry3.config(state=tkinter.NORMAL)
           self.label12.config(state=tkinter.NORMAL)
           self.entry12.config(state=tkinter.NORMAL)
        self.adjust_parameters()

    def disable_geometry(self):
        self.label32.config(state=tkinter.DISABLED)
        self.entry32.config(state=tkinter.DISABLED)
        self.label33.config(state=tkinter.DISABLED)
        self.entry33.config(state=tkinter.DISABLED)
        self.label34.config(state=tkinter.DISABLED)
        self.entry34.config(state=tkinter.DISABLED)
        self.label35.config(state=tkinter.DISABLED)
        self.entry35.config(state=tkinter.DISABLED)
        self.label36.config(state=tkinter.DISABLED)
        self.entry36.config(state=tkinter.DISABLED)
        self.label37.config(state=tkinter.DISABLED)
        self.entry37.config(state=tkinter.DISABLED)
        self.label38.config(state=tkinter.DISABLED)
        self.entry38.config(state=tkinter.DISABLED)
        self.label39.config(state=tkinter.DISABLED)
        self.entry39.config(state=tkinter.DISABLED)
        self.label40.config(state=tkinter.DISABLED)
        self.engine.config(state=tkinter.DISABLED)
        self.label41.config(state=tkinter.DISABLED)        
        self.entry41.config(state=tkinter.DISABLED)
        self.label42.config(state=tkinter.DISABLED)        
        self.entry42.config(state=tkinter.DISABLED)
        self.label43.config(state=tkinter.DISABLED)        
        self.entry43.config(state=tkinter.DISABLED)
        self.label44.config(state=tkinter.DISABLED)        
        self.entry44.config(state=tkinter.DISABLED)
        self.label45.config(state=tkinter.DISABLED)      
        self.entry45.config(state=tkinter.DISABLED)
        self.cg_check.config(state=tkinter.DISABLED)  
        self.label46.config(state=tkinter.DISABLED)      
        self.entry46.config(state=tkinter.DISABLED)
        self.label47.config(state=tkinter.DISABLED)        
        self.entry47.config(state=tkinter.DISABLED)
        self.label48.config(state=tkinter.DISABLED)        
        self.entry48.config(state=tkinter.DISABLED)
        self.label49.config(state=tkinter.DISABLED)        
        self.entry49.config(state=tkinter.DISABLED)
        self.label50.config(state=tkinter.DISABLED)        
        self.entry50.config(state=tkinter.DISABLED)
        self.label51.config(state=tkinter.DISABLED)        
        self.entry51.config(state=tkinter.DISABLED)
        self.label52.config(state=tkinter.DISABLED)        
        self.entry52.config(state=tkinter.DISABLED)

    def gsolver_change(self,*args):
        self.disable_geometry()
        if self.solver_type.get() == 'Minimal':
           self.label32.config(state=tkinter.NORMAL)
           self.entry32.config(state=tkinter.NORMAL)
        elif self.solver_type.get() == 'Evolutionary':
           self.label33.config(state=tkinter.NORMAL)
           self.entry33.config(state=tkinter.NORMAL)
           self.label34.config(state=tkinter.NORMAL)
           self.entry34.config(state=tkinter.NORMAL)
           self.label35.config(state=tkinter.NORMAL)
           self.entry35.config(state=tkinter.NORMAL)
        elif self.solver_type.get() == 'Annealing':
           self.label36.config(state=tkinter.NORMAL)
           self.entry36.config(state=tkinter.NORMAL)
           self.label37.config(state=tkinter.NORMAL)
           self.entry37.config(state=tkinter.NORMAL)
           self.label38.config(state=tkinter.NORMAL)
           self.entry38.config(state=tkinter.NORMAL)
           self.label39.config(state=tkinter.NORMAL)
           self.entry39.config(state=tkinter.NORMAL)
        elif self.solver_type.get() == 'Mechanical':
           self.label40.config(state=tkinter.NORMAL)
           self.engine.config(state=tkinter.NORMAL)
           self.cg_check.config(state=tkinter.NORMAL)
           self.label41.config(state=tkinter.NORMAL)
           self.entry41.config(state=tkinter.NORMAL)
           self.label42.config(state=tkinter.NORMAL)
           self.entry42.config(state=tkinter.NORMAL)
           self.label43.config(state=tkinter.NORMAL)
           self.entry43.config(state=tkinter.NORMAL)
           self.label44.config(state=tkinter.NORMAL)
           self.entry44.config(state=tkinter.NORMAL)
           self.label45.config(state=tkinter.NORMAL)
           self.entry45.config(state=tkinter.NORMAL)
           self.label46.config(state=tkinter.NORMAL)
           self.entry46.config(state=tkinter.NORMAL)
           self.label47.config(state=tkinter.NORMAL)
           self.entry47.config(state=tkinter.NORMAL)
           self.label48.config(state=tkinter.NORMAL)
           self.entry48.config(state=tkinter.NORMAL)
        else:
           self.label49.config(state=tkinter.NORMAL)
           self.entry49.config(state=tkinter.NORMAL)
           self.label50.config(state=tkinter.NORMAL)
           self.entry50.config(state=tkinter.NORMAL)
           self.label51.config(state=tkinter.NORMAL)
           self.entry51.config(state=tkinter.NORMAL)
           self.label52.config(state=tkinter.NORMAL)
           self.entry52.config(state=tkinter.NORMAL)

    def load_parameters(self):
        if (self.parameter_filename.get() == ""):
            filename = askopenfilename(filetypes=(("XML File", "*.xml"),("All Files","*.*")))
            if filename:
            	self.parameter_filename.set(os.path.basename(filename))
            else:
                return
            self.read_parameters()

    def save_parameters(self):
        if (self.parameter_filename.get() == ""):
            filename = tkinter.filedialog.SaveFileDialog(root).go("*.xml")
            if filename:
                self.parameter_filename.set(os.path.basename(filename))
            else:
                return
        self.write_parameters()

    def clear_parameters(self):
        self.sheet_dynamics.set(False)
        self.superposable.set(True)
        self.compressible.set(True)
        self.permutable.set(True)
        self.perturbg.set(False)
        self.perturbt.set(True)
        self.perturbe.set(True)
        self.dim_uniformity.set(True)
        self.relational.set(False)
        self.high_memory_footprint.set(True)
        self.cg_refinement.set(True)

        self.initial_state.set('Cartesian')
        self.hyphansis.set('Dynamic')
        self.homology_method.set('Native')
        self.homology_field.set('INT')
        self.signature.set('Lorentzian')
        self.solver_type.set('Minimal')
        self.int_engine.set('RK4')

        self.initial_events.set(1296)
        self.max_iterations.set(25)
        self.initial_dimension.set(4)
        self.initial_sheets.set(1)
        self.chkpt_frequency.set(5)
        self.background_dim.set(4)
        self.random_seed.set(0)
        self.solver_iterations.set(100)
        self.max_generations.set(100)
        self.pool_size.set(50)
        self.max_jousts.set(20)
        self.thermal_sweeps.set(1000)
        self.annealing_steps.set(500)
        self.max_int_steps.set(10000)
        self.max_cg_steps.set(20)
        self.max_ls_steps.set(10)

        self.edge_probability.set(0.15)
        self.abnormality_threshold.set(0.1)
        self.superposition_threshold.set(0.05)
        self.parity_probability.set(0.05)
        self.geometry_threshold.set(0.00005)
        self.thermal_variance.set(0.5)
        self.thermalization_criterion.set(0.001)
        self.step_size.set(0.05)
        self.damping.set(0.85)
        self.spring.set(-1.5)
        self.repulsion.set(1.0)
        self.edge_flexibility.set(2.0)
        self.reflection.set(0.9)
        self.contraction.set(2.2)
        self.expansion.set(0.75)
        self.shrinkage.set(0.45)
        
    def convert_boolean(self,bvalue):
        if bvalue == 0:
            return "No"
        else:
            return "Yes"

    def convert_string(self,val):
        val.upper()
        if val == "YES":
            return 1
        else:
            return 0

    def adjust_parameters(self):
        if self.initial_state.get() == 'Cartesian':
            self.initial_events.set(1296)
            self.initial_dimension.set(4)
        if self.initial_state.get() == 'Singleton':
            self.initial_events.set(1)
            self.initial_dimension.set(0)
        elif self.initial_state.get() == 'Monoplex':
            self.initial_events.set(5)
            self.initial_dimension.set(4)
        elif self.initial_state.get() == 'Random':
            self.initial_events.set(1000)
            self.initial_dimension.set(4)

    def initial_event_change(self,*args):
        if self.initial_state.get() == 'Cartesian':
            n = self.initial_events.get()
            d = self.initial_dimension.get()
            self.initial_events.set(int(math.pow(n,1.0/double(d))))

    def initial_dimension_change(self,*args):
        if self.initial_state.get() == 'Cartesian':
            n = self.initial_events.get()
            d = self.initial_dimension.get()
            self.initial_events.set(int(math.pow(n,1.0/double(d))))
        elif self.initial_state.get() == 'Monoplex':
            self.initial_events.set(1+self.initial_dimension.get())

    def read_parameters(self):
        parameter_filename = self.parameter_filename.get()

    def write_parameters(self):
        parameter_filename = self.parameter_filename.get()
        
root = tkinter.Tk()
gui = euplecton(root)
root.mainloop()










