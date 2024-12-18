import numpy as np
import matplotlib.pyplot as plt

"""
O código abaixo apenas plota a trajetória determinada por diferentes métodos numéricos sujeitas a uma pequena perturbação de 0.01 nas 
condiçoes iniciais, descritas em estado_inicial.
"""

# --- Adicionando perturbações nas condições iniciais ---
perturbacao = np.array([0.01, 0.01, 0.01])

# --- Parâmetros globais ---
sigma = 10.0
rho = 28.0
beta = 8 / 3
t0, tf = 0, 50  # Tempo inicial e final
h = 0.001       # Passo de tempo
num_steps = int((tf - t0) / h)  # Número total de passos
eta = np.sqrt(beta * (rho - 1))  # Valor de eta
estado_inicial = np.array([rho - 1, eta, eta - 3])  # Condições iniciais

# --- Método de Euler ---
def metodo_euler(estado, h, num_steps):
    x, y, z = estado
    trajetoria = [estado]
    for _ in range(num_steps):
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z
        x += h * dx
        y += h * dy
        z += h * dz
        trajetoria.append([x, y, z])
    return np.array(trajetoria)

# --- Método de Heun ---
def metodo_heun(estado, h, num_steps):
    x, y, z = estado
    trajetoria = [estado]
    for _ in range(num_steps):
        # Predição inicial (Euler)
        dx_euler = sigma * (y - x)
        dy_euler = x * (rho - z) - y
        dz_euler = x * y - beta * z
        x_pred = x + h * dx_euler
        y_pred = y + h * dy_euler
        z_pred = z + h * dz_euler

        # Correção (Heun)
        dx_corr = sigma * (y_pred - x_pred)
        dy_corr = x_pred * (rho - z_pred) - y_pred
        dz_corr = x_pred * y_pred - beta * z_pred
        x += h / 2 * (dx_euler + dx_corr)
        y += h / 2 * (dy_euler + dy_corr)
        z += h / 2 * (dz_euler + dz_corr)
        trajetoria.append([x, y, z])
    return np.array(trajetoria)

# --- Método RK4 ---
def sistema_lorenz(estado, sigma, rho, beta):
    x, y, z = estado
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return np.array([dx, dy, dz])

def runge_kutta_4(estado, h, num_steps):
    trajetoria = [estado]
    for _ in range(num_steps):
        k1 = h * sistema_lorenz(estado, sigma, rho, beta)
        k2 = h * sistema_lorenz(estado + 0.5 * k1, sigma, rho, beta)
        k3 = h * sistema_lorenz(estado + 0.5 * k2, sigma, rho, beta)
        k4 = h * sistema_lorenz(estado + k3, sigma, rho, beta)
        estado = estado + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        trajetoria.append(estado)
    return np.array(trajetoria)

# --- Programa principal ---
trajetoria_euler = metodo_euler(estado_inicial.copy(), h, num_steps)
trajetoria_heun = metodo_heun(estado_inicial.copy(), h, num_steps)
trajetoria_rk4 = runge_kutta_4(estado_inicial.copy(), h, num_steps)

# Novas condições iniciais com perturbações
estado_inicial_perturbado = estado_inicial + perturbacao

# Calculando as trajetórias para as condições perturbadas
trajetoria_euler_perturbada = metodo_euler(estado_inicial_perturbado.copy(), h, num_steps)
trajetoria_heun_perturbada = metodo_heun(estado_inicial_perturbado.copy(), h, num_steps)
trajetoria_rk4_perturbada = runge_kutta_4(estado_inicial_perturbado.copy(), h, num_steps)

# --- Calculando a diferença nas trajetórias ---
diferenca_euler = np.linalg.norm(trajetoria_euler - trajetoria_euler_perturbada, axis=1)
diferenca_heun = np.linalg.norm(trajetoria_heun - trajetoria_heun_perturbada, axis=1)
diferenca_rk4 = np.linalg.norm(trajetoria_rk4 - trajetoria_rk4_perturbada, axis=1)

# --- Visualização ---
fig = plt.figure(figsize=(17, 8))
plt.plot(np.linspace(t0, tf, num_steps + 1), diferenca_euler, label='Diferença Euler', color='red')
plt.plot(np.linspace(t0, tf, num_steps + 1), diferenca_heun, label='Diferença Heun', color='blue')
plt.plot(np.linspace(t0, tf, num_steps + 1), diferenca_rk4, label='Diferença RK4', color='green')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.title('Diferença nas Trajetórias', fontsize=22)
plt.xlabel('Tempo', fontsize=20)
plt.ylabel('Norma da Diferença', fontsize=20)
plt.legend(fontsize=20)
plt.tight_layout()
plt.show()
