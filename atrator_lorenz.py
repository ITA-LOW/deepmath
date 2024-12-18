import numpy as np
import matplotlib.pyplot as plt

"""
O código abaixo apenas plota a trajetória determinada por diferentes métodos numéricos.
As análises solicitadas na entrevista estão nos demais códigos. Aumente o valor da variável 'tf' na linha 14 
para obter uma representação mais próxima da borboleta que dá nome ao efeito caótico gerado pelo atrator de Lorenz.
"""

# --- Parâmetros globais ---
sigma = 10.0
rho = 28.0
beta = 8 / 3
t0, tf = 0, 10  # Tempo inicial e final
h = 0.001       # Passo de tempo
num_passos = int((tf - t0) / h)  # Número total de passos
eta = np.sqrt(beta * (rho - 1))  # Valor de eta
estado_inicial = np.array([rho - 1, eta, eta - 3])  # Condições iniciais


# --- Método de Euler ---
def metodo_euler(estado, h, num_passos):
    x, y, z = estado
    trajetoria = [estado]
    for _ in range(num_passos):
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z
        x += h * dx
        y += h * dy
        z += h * dz
        trajetoria.append([x, y, z])
    return np.array(trajetoria)


# --- Método de Heun ---
def metodo_heun(estado, h, num_passos):
    x, y, z = estado
    trajetoria = [estado]
    for _ in range(num_passos):
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


def runge_kutta_4(estado, h, num_passos):
    trajetoria = [estado]
    for _ in range(num_passos):
        k1 = h * sistema_lorenz(estado, sigma, rho, beta)
        k2 = h * sistema_lorenz(estado + 0.5 * k1, sigma, rho, beta)
        k3 = h * sistema_lorenz(estado + 0.5 * k2, sigma, rho, beta)
        k4 = h * sistema_lorenz(estado + k3, sigma, rho, beta)
        estado = estado + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        trajetoria.append(estado)
    return np.array(trajetoria)


# --- Programa principal ---
trajetoria_euler = metodo_euler(estado_inicial.copy(), h, num_passos)
trajetoria_heun = metodo_heun(estado_inicial.copy(), h, num_passos)
trajetoria_rk4 = runge_kutta_4(estado_inicial.copy(), h, num_passos)


# --- Visualização ---
fig = plt.figure(figsize=(17, 8))

# Plot do Método de Euler
ax1 = fig.add_subplot(131, projection='3d')
ax1.plot(trajetoria_euler[:, 0], trajetoria_euler[:, 1], trajetoria_euler[:, 2], lw=1.2, color="red")
ax1.set_title("Método de Euler", fontsize=22)
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.set_zlabel("Z")

# Plot do Método de Heun
ax2 = fig.add_subplot(132, projection='3d')
ax2.plot(trajetoria_heun[:, 0], trajetoria_heun[:, 1], trajetoria_heun[:, 2], lw=1.2, color="blue")
ax2.set_title("Método de Heun", fontsize=22)
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.set_zlabel("Z")

# Plot do Método RK4
ax3 = fig.add_subplot(133, projection='3d')
ax3.plot(trajetoria_rk4[:, 0], trajetoria_rk4[:, 1], trajetoria_rk4[:, 2], lw=1.2, color="green")
ax3.set_title("Método RK4", fontsize=22)
ax3.set_xlabel("X")
ax3.set_ylabel("Y")
ax3.set_zlabel("Z")

plt.tight_layout()
plt.show()
