#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#define FPS 60
#define DT (1.0f / FPS)
#define N_STEPS 100
#define SDT (DT / N_STEPS)
#define N_BEADS 5

typedef float f32;
typedef int32_t i32;

typedef struct vec2 {
  f32 x;
  f32 y;
} vec2;

static vec2 gravity = {0.0F, -10.0F};

typedef struct bead {
  f32 radius;
  f32 mass;
  vec2 pos;
  vec2 prev_pos;
  vec2 vel;
} bead;

typedef struct wire {
  vec2 pos;
  f32 radius;
} wire;

static bead bds[N_BEADS]; 
static wire wr = {{0.0f, 0.0f}, 0.8f};

static vec2 vec2_sub(vec2 a, vec2 b) {
  return (vec2) {a.x - b.x, a.y - b.y};
}

static vec2 vec2_adds(vec2 a, vec2 b, float s) {
  return (vec2) {a.x + b.x * s, a.y + b.y * s};
}

static vec2 vec2_subs(vec2 a, vec2 b, float s) {
  return (vec2) {a.x - b.x * s, a.y - b.y * s};
}

static f32 vec2_dot(vec2 a, vec2 b) {
  return a.x * b.x + a.y * b.y;
}

static vec2 vec2_divs(vec2 a, float b) {
  return (vec2) {a.x / b, a.y / b};
}

static f32 vec2_len(vec2 a) {
  return sqrtf(vec2_dot(a, a));
}

static void start_step(bead *bd) {
  bd->vel = vec2_adds(bd->vel, gravity, SDT);
  bd->prev_pos = bd->pos;
  bd->pos = vec2_adds(bd->pos, bd->vel, SDT);
}

static void end_step(bead *bd) {
  bd->vel = vec2_sub(bd->pos, bd->prev_pos);
  bd->vel = vec2_divs(bd->vel, SDT);
}

static void bead_col(bead *a, bead *b) {
  vec2 dir = vec2_sub(b->pos, a->pos);
  f32 d = vec2_len(dir);
  if (d == 0.0f || d > a->radius + b->radius)
    return;
  dir = vec2_divs(dir, d);
  f32 corr = (a->radius + b->radius - d) / 2.0f;
  a->pos = vec2_subs(a->pos, dir, corr);
  b->pos = vec2_adds(b->pos, dir, corr);
  f32 v0a = vec2_dot(a->vel, dir); 
  f32 v0b = vec2_dot(b->vel, dir); 
  f32 ma = a->mass;
  f32 mb = b->mass;
  f32 mt = ma + mb;
  f32 vc = ma * v0a + mb * v0b;
  f32 v1a = (vc - mb * (v0a - v0b)) / mt;
  f32 v1b = (vc - ma * (v0b - v0a)) / mt; 
  a->vel = vec2_adds(a->vel, dir, v1a - v0a);
  b->vel = vec2_adds(b->vel, dir, v1b - v0b);
}

static void keep_on_wire(bead *bd, wire *wr) {
  vec2 dir = vec2_sub(bd->pos, wr->pos);
  f32 len = vec2_len(dir); 
  if (len == 0.0f)
    return;
  dir = vec2_divs(dir, len);
  f32 lambda = wr->radius - len;
  bd->pos = vec2_adds(bd->pos, dir, lambda);
}

static void print_sim(int f) {
  for (i32 i = 0; i < N_BEADS; i++) {
    bead *bd = bds + i;
    printf("%d,%d,%f,%f,%f\n", f, 0, bd->pos.x, bd->pos.y, bd->radius); 
  }
  printf("%d,%d,%f,%f,%f\n", f, 1, wr.pos.x, wr.pos.y, wr.radius); 
}

int main(int argc, char **argv) {
  srand48(time(NULL));
  switch (argc) {
  case 0:
  case 1:
    break;
  case 2:
    FILE *tmp;
    tmp = freopen(argv[1], "wb", stdout);
    if (!tmp) {
      perror("freopen");
      exit(1);
    }
    stdout = tmp;
    break;
  default:
    fputs("too many args\n", stderr);
    exit(1);
  }
  f32 r = 0.1f;
  f32 rot = 0.0f;
  for (i32 i = 0; i < N_BEADS; i++) {
    bead *bd = bds + i;
    bd->radius = r;
    bd->mass = (float) M_PI * r * r; 
    bd->pos.x = wr.pos.x + wr.radius * cosf(rot);
    bd->pos.y = wr.pos.y + wr.radius * sinf(rot);
    rot += (float) M_PI / N_BEADS;
    r = 0.05f + drand48() * 0.1f;
  }
  printf("f,t,x,y,r\n");
  i32 f;
  for (f = 0; f < FPS * 10; f++) {
    print_sim(f);
    for (i32 s = 0; s < N_STEPS; s++) {
      i32 i, j;
      for (i = 0; i < N_BEADS; i++)
        start_step(bds + i);
      for (i = 0; i < N_BEADS; i++)
        keep_on_wire(bds + i, &wr);
      for (i = 0; i < N_BEADS; i++)
        end_step(bds + i);
      for (i = 0; i < N_BEADS; i++) {
        for (j = 0; j < i; j++)
          bead_col(bds + i, bds + j);
      }
    }
  }
  print_sim(f);
}
