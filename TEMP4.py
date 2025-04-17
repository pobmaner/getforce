import tkinter as tk
from tkinter import messagebox
import math
import numpy as np
from sympy import symbols, solve
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
e=0.000000000000002
def xtrany(z, f):
    return z * f

def xtranz(y, f):
    return y * f

def ytranz(x, f):
    return x * f

def ytranx(z, f):
    return z * f

def ztranx(y, f):
    return y * f

def ztrany(x, f):
    return x * f

def dep2(x, y, z):
    return x * x + y * y + z * z

def depsom2(a, b, c, x, y, z):
    return (a - x) * (a - x) + (b - y) * (b - y) + (c - z) * (c - z)
def dep2(x, y, z):
    return x * x + y * y + z * z

def depsom2(a, b, c, x, y, z):
    return (a - x) * (a - x) + (b - y) * (b - y) + (c - z) * (c - z)

def jiacosta(fx, fy, fz, mx, my, mz):
    costa = (dep2(fx, fy, fz) + dep2(mx, my, mz) - depsom2(fx, fy, fz, mx, my, mz)) / (2 * math.sqrt(dep2(fx, fy, fz)) * math.sqrt(dep2(mx, my, mz)))
    if costa > 0:
        return 1
    else:
        return 0


root = tk.Tk()
root.title("力系简化程序 2024.10")


tk.Label(root, text="力的数量:").grid(row=0, column=0)
n_entry = tk.Entry(root)
n_entry.grid(row=0, column=1)

tk.Label(root, text="该力的 x 分量:(米)").grid(row=1, column=0)
x1_entry = tk.Entry(root)
x1_entry.grid(row=1, column=1)

tk.Label(root, text="该力的 x 分量的大小:(N)").grid(row=2, column=0)
fx_entry = tk.Entry(root)
fx_entry.grid(row=2, column=1)

tk.Label(root, text="该力的 y 分量:(米)").grid(row=3, column=0)
y1_entry = tk.Entry(root)
y1_entry.grid(row=3, column=1)

tk.Label(root, text="该力的 y 分量的大小:(N)").grid(row=4, column=0)
fy_entry = tk.Entry(root)
fy_entry.grid(row=4, column=1)

tk.Label(root, text="该力的 z 分量:(米)").grid(row=5, column=0)
z1_entry = tk.Entry(root)
z1_entry.grid(row=5, column=1)

tk.Label(root, text="该力的 z 分量的大小:(N)").grid(row=6, column=0)
fz_entry = tk.Entry(root)
fz_entry.grid(row=6, column=1)


fx, fy, fz = 0, 0, 0
mx, my, mz = 0, 0, 0
cue = 0


result_text = tk.Text(root, height=10, width=50)
result_text.grid(row=8, column=0, columnspan=2)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=9, column=0, columnspan=2)


def calculate():
    global fx, fy, fz, mx, my, mz, cue
    try:
        n = int(n_entry.get())
        x1 = float(x1_entry.get())
        fx1 = float(fx_entry.get())
        y1 = float(y1_entry.get())
        fy1 = float(fy_entry.get())
        z1 = float(z1_entry.get())
        fz1 = float(fz_entry.get())

        fx += fx1
        fy += fy1
        fz += fz1
        mx +=  ztranx(y1, fz1)-ytranx(z1, fy1)
        my += xtrany(z1, fx1) - ztrany(x1, fz1)
        mz += ytranz(x1, fy1)-xtranz(y1, fx1)  


        cue += 1

        if cue >= n:
            show_result()
        else:
            preview_result()
            clear_entries()
    except Exception as e:
        messagebox.showerror("错误", str(e))

def clear_entries():
    x1_entry.delete(0, tk.END)
    fx_entry.delete(0, tk.END)
    y1_entry.delete(0, tk.END)
    fy_entry.delete(0, tk.END)
    z1_entry.delete(0, tk.END)
    fz_entry.delete(0, tk.END)

def preview_result():
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, f"当前主矢：Fx={fx}, Fy={fy}, Fz={fz}\n")
    result_text.insert(tk.END, f"当前主矩：Mx={mx}, My={my}, Mz={mz}\n")
    update_plot()

def update_plot():
    ax.clear()
    ax.quiver(0, 0, 0, fx, fy, fz, color='r', label='Principal vector')
    ax.quiver(0, 0, 0, mx, my, mz, color='b', label='Moment')
    ax.set_xlim([-100,100])
    ax.set_ylim([-100,100])
    ax.set_zlim([-100,100])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    canvas.draw()

def show_result():
    dff = dep2(fx, fy, fz)
    dmm = dep2(mx, my, mz)

    if dff == 0 or dmm == 0:
        result = "主矢：Fx={}, Fy={}, Fz={}\n主矩：Mx={}, My={}, Mz={}\n".format(fx, fy, fz, mx, my, mz)
        if dep2(fx, fy, fz) + dep2(mx, my, mz) == 0:
            result += "简化为平衡力系"
        elif dep2(fx, fy, fz) == 0:
            result += "简化为力偶系"
        else:
            result += "简化为合力"
    else:
        costa = (dep2(fx, fy, fz) + dep2(mx, my, mz) - depsom2(fx, fy, fz, mx, my, mz)) / (2 * math.sqrt(dep2(fx, fy, fz)) * math.sqrt(dep2(mx, my, mz)))
        if abs(costa-1)<e:
            costa=1
            sinta=0
            dtodd=0
        elif abs(costa+1)<e:
            costa=-1
            sinta=0
            dtodd=0
        else:
            if abs(costa)<e:
                costa=0
                sinta=1
            else:
                sinta = math.sqrt(1 - costa * costa)
            cha1 = []
            ff = [fx, fy, fz]
            mm = [mx, my, mz]
            cha1 = np.cross(ff, mm)
        
            a, b, c = symbols('a b c')

            costesta = (dep2(fx, fy, fz) + dep2(mx, my, mz) - depsom2(fx, fy, fz, mx, my, mz)) / (2 * math.sqrt(dep2(fx, fy, fz)) * math.sqrt(dep2(mx, my, mz)))
            mm2two = solve([cha1[0] * a + cha1[1] * b + cha1[2] * c, a * a + b * b + c * c - dep2(mx, my, mz) * sinta * sinta, fx * a + fy * b + fz * c], [a, b, c])
            if sinta==1:
                mm2=mm
            elif jiacosta(mm2two[0][0], mm2two[0][1], mm2two[0][2], mx, my, mz) == 1:
                mm2 = mm2two[0]
            else:
                mm2 = mm2two[1]
        
            dd = np.cross(ff, mm2)
            dmm2 = dep2(mm2[0], mm2[1], mm2[2])
            dtodd = dep2(dd[0], dd[1], dd[2])
            d = math.sqrt(dmm2) / math.sqrt(dff)
            
            k=np.dot([ff[0],ff[1],ff[2]],[mm[0],mm[1],mm[2]])/dff
            lx=[1,1,1]
            for i in range(3):
                lx[i]=k*ff[i]
            
            for i in range(len(dd)):
                dd[i] /=dep2(ff[0] , ff[1] , ff[2])
            dtodd = dep2(dd[0], dd[1], dd[2])
        result = "主矢：Fx={}, Fy={}, Fz={}\n主矩：Mx={}, My={}, Mz={}\n".format(fx, fy, fz, mx, my, mz)
            
        if abs(dtodd)<e:
            result += "简化为过简化点的力螺旋，力螺旋的矩为主矩"
        else:
            if costa == 0:
                result += "简化为不过简化点的合力，作用点的位置向量为({},{},{})".format(dd[0], dd[1], dd[2])
            else:
                result += "简化为不过简化点的力螺旋，作用点的位置向量为({},{},{})".format(dd[0], dd[1], dd[2])
                result += "，力螺旋的矩向量为({},{},{})".format(lx[0], lx[1], lx[2])
    
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, result)
    update_plot()


tk.Button(root, text="提交当前力", command=calculate).grid(row=7, column=0, columnspan=2)


root.mainloop()