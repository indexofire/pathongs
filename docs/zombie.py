import itertools, sys, time, random, math, pygame
from pygame.locals import *
from MyLibrary import *

def calc_velocity(direction, vel=1.0):
    velocity = point(0, 0)
    if direction == 0:
        velocity.y = -vel
    elif direction == 2:
        velocity.x = vel
    elif direction == 4:
        velocity.y = vel
    elif direction == 6:
        velocity.x = -vel
    return velocity

def reverse_direction(sprite):
    if sprite.direction == 0:
        sprite.direction = 4
    elif sprite.direction == 2:
        sprite.direction = 6
    elif sprite.direction == 4:
        sprite.direction = 0
    elif sprite.direction == 6:
        sprite.direction = 2

pygame.init()
screen = pygame.display.set_mode((800,600))
pygame.display.set_caption("collision demo yuguo")
font = pygame.font.Font(None, 36)
timer = pygame.time.Clock()

player_group = pygame.sprite.Group()
zombie_group = pygame.sprite.Group()
health_group = pygame.sprite.Group()

player = MySprite()
player.load("farmer walk.png", 96, 96, 8)
player.position = 80, 80
player.direction = 4
player_group.add(player)
