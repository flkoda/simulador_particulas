import pygame, random
import numpy as np
import pygame.gfxdraw

class Ambiente():
    def __init__(self, DIM, dt):
        self.DIM = DIM
        self.dt = dt
        self.particles = []

    def update(self):
        for p1 in self.particles:
            p1.stateUpdate()
            self.bounce(p1)
            for p2 in self.particles:
                if p1 != p2:
                    self.elasticCollision(p1, p2)

    def addParticle(self, p):
        self.particles.append(p)

    def bounce(self, p):
        for i in range(2):
            if p.X[0][i] > self.DIM[i] - p.radius:
                p.X[0][i] = self.DIM[i] - p.radius
                p.V[0][i] *= -1
            elif p.X[0][i] < p.radius:
                p.X[0][i] = p.radius
                p.V[0][i] *= -1

    def elasticCollision(self, p1, p2):
        dX = p2.X - p1.X
        dist = np.linalg.norm(dX)

        if dist < p1.radius + p2.radius:
            overlap = 0.5 * (p1.radius + p2.radius - dist)
            p1.X -= dX * overlap / dist
            p2.X += dX * overlap / dist

            n = dX / dist
            dv = p1.V - p2.V

           
            impulse_scalar = (2 * np.dot(dv[0], n[0])) / (p1.mass + p2.mass)  # Extrai os vetores 1D

            # Aplica o impulso como vetores
            impulse = impulse_scalar * n

            p1.V -= p2.mass * impulse
            p2.V += p1.mass * impulse


class Particula():
    def __init__(self, env, X, V, radius, mass):
        self.env = env
        self.X = X
        self.V = V
        self.radius = radius
        self.mass = mass
        self.colour = (0, int(random.uniform(0,255)), 0)

    def stateUpdate(self):
        self.X += self.V * self.env.dt


# --- Inicialização ---
DIM = np.array([700, 700])
dt = 0.01
env = Ambiente(DIM, dt)

pygame.init()
screen = pygame.display.set_mode((DIM[0], DIM[1]))
pygame.display.set_caption('Simulador ')
number_of_particles = 50

for _ in range(number_of_particles):
    radius = np.random.randint(10, 20)
    mass = radius**3
    X = np.random.rand(1, 2) * (DIM - radius) + radius
    V = (np.random.rand(1, 2) - 0.5) * 100

    # Certifique-se que X e V são arrays NumPy com tipo float
    X = X.astype(float)
    V = V.astype(float)
    particle = Particula(env, X, V, radius, mass)
    env.addParticle(particle)

def display(env):
    for p in env.particles:
        pygame.gfxdraw.filled_circle(screen, int(p.X[0][0]), int(p.X[0][1]), p.radius, p.colour)


# --- Loop principal ---
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill([255, 255, 255])
    env.update()
    display(env)
    pygame.display.flip()

pygame.quit()

