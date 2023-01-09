# A script that attempts ray tracing. Sites that have helped construct this are: 
# https://medium.com/swlh/ray-tracing-from-scratch-in-python-41670e6a96f9
# https://raytracing.github.io/books/RayTracingInOneWeekend.html#surfacenormalsandmultipleobjects

# Harun Gürhan hgurhan@kth.se
# Last updated 2022-12-12

#Libraries ---- 
import math
import numpy as np
from tkinter import *
from tkinter import ttk
import numpy as np
import sys

image_width = 480
image_height = 270


def ray_tracing_algorithm(camera=np.array([0, 0, 0])):
    """Creates and sends out light rays from a starting point, the eye/camera
    Calls the nearest object function to see if they intersect any object in space
    If they do intersect, calls the reflection function
    IN: Point in coordinate space
    OUT: Multiple light rays."""

    global frame
    #Sets up the screen
    tally = 0
    upper_left=np.array([-(16.0/9.0)*1, 1, -1])
    #cli_color = None
    for x in range(image_width):
        #print(cli_color)
        #cli_color = ""
        i = x/image_width
        for y in range(image_height):
            j = y/image_height
            frame.img.put("black", (x, y))
            #calculates which pixel to calculate color of
            current_pixel = upper_left + i*np.array([(16/9)*2, 0, 0]) - j*np.array([0, 2, 0])

            d_vec = normalize_vector(current_pixel - camera)
            tally += 1
            light_ray = Light(start_point=camera, vec=d_vec)
            bounces = 0
            closest = nearest_intersected_obj(light_ray)
            c_p = 0
            reflectionco = 1
            if closest != None:
                closest_object, min_distance = closest[0], closest[1]
                intersection_point = light_ray.start_point + min_distance*light_ray.vec
                normal = normalize_vector(intersection_point - closest_object.centerpoint)
                #Checks if shadow
                L = normalize_vector(light_source.point - intersection_point + 0.00001*normal)
                ray_to_light = Light(start_point=intersection_point + 0.00001*normal, vec=L)
                next_ = nearest_intersected_obj(ray_to_light) 
                #if not shadow, calculate color of current pixel with reflections
                if next_ == None or next_ != None and np.linalg.norm(L) < next_[1]:
                    #Checks if ray comes from outside the sphere or inside
                    front_face = ray_from_outside_sphere(normal, light_ray.vec)
                    if front_face == True:
                        while bounces < 5:
                            reflection_ray = reflection(intersection_point, closest_object, d_vec)
                            next_ = nearest_intersected_obj(reflection_ray)
                            
                            if bounces == 0:
                                c_p += calculate_color(camera=camera, intersection_point=intersection_point, obj=closest_object)
                                reflectionco *= closest_object.reflectionco
                            elif bounces > 0:
                                c_p += reflectionco*calculate_color(camera=camera, intersection_point=intersection_point, obj=closest_object)                    
                                reflectionco *= closest_object.reflectionco

                            if next_ == None:
                                bounces = 6
                    
                            elif next_ != None:
                                closest = next_
                                closest_object, min_distance = closest[0], closest[1]
                                light_ray = reflection_ray
                                intersection_point = light_ray.start_point + min_distance*light_ray.vec

                                bounces += 1          
                    elif front_face == False:
                        continue
                elif next_ != None:
                    continue
            
            #for cli print
            """if type(c_p) == type(0):
                    cli_color += " M" 
            """

            if type(c_p) != type(0):
                for s in range(0, 3):
                    if c_p[s] > 1:
                        c_p[s] = 1
                    if c_p[s] < 0:
                        c_p[s] = 0
                
                #for cli print
                """if 0 < c_p[1] <= 0.3:
                    cli_color += " *"
                if 0.3 < c_p[1] <= 0.5:
                    cli_color += " +"
                if 0.5 < c_p[1] <= 0.7:
                    cli_color += " -"
                if 0.7 < c_p[1] <= 0.9:
                    cli_color += " ."
                if 0.9 < c_p[1] <= 1:
                    cli_color += " l"   """

                
                # translates from rgb in range [0, 1] to rgb in range [0, 255] to hex color for tkinter
                
                c_p = (abs(int(255*c_p[0])), abs(int(255*c_p[1])), abs(int(255*c_p[2])))
                c_p = rgbtohex(c_p[0], c_p[1], c_p[2])

                #progress counter https://stackoverflow.com/questions/6169217/replace-console-output-in-python
                progress = 100.0*(tally/(image_height*image_width))
                sys.stdout.write("\rprogress: %d/100" % progress)
                sys.stdout.flush
                frame.img.put(c_p, (x, y))
                    
def normalize_vector(vec):
    """"""
    return np.divide(vec, np.linalg.norm(vec))

def nearest_intersected_obj(light_ray):
    """Goes through objects in scene to see which ones the ray intersects. 
    Counts only the closest one. If the ray does not intersect anything returns None.
    """
    intersected_obj_list = []
    for obj in range(len(objects)):
        intersection_ = intersection(light_ray, objects[obj])
        if intersection_ != None:
            intersected_obj_list.append(intersection_)
            intersected_obj_list.append(objects[obj])
    if len(intersected_obj_list) > 0:
        min_distance, closest_object = intersected_obj_list[0], intersected_obj_list[1]
        for l in range(0, len(intersected_obj_list), 2):
            for k in range(0, len(intersected_obj_list), 2):
                if intersected_obj_list[l] < intersected_obj_list[k] and intersected_obj_list[l] < min_distance:
                    min_distance, closest_object  = intersected_obj_list[l], intersected_obj_list[l + 1]
        return (closest_object, min_distance)
    else:
        return None

def intersection(light_ray, sphere):
    """Solves for t in the equation abs(t*direction_vector + origin) = radius of sphere
    in order to see if there is a point in space where the vector intersects the sphere
    Returns t if there is or None if there is not"""

    a = (np.linalg.norm(light_ray.vec))**2
    b = 2*np.dot(light_ray.vec, (light_ray.start_point - sphere.centerpoint))
    c = (np.linalg.norm(light_ray.start_point - sphere.centerpoint))**2 - sphere.radius**2
    discriminant = b**2 - 4*a*c
    if discriminant <= 0:
        return None
    t1 = ((-b + math.sqrt(discriminant))/(2*a))
    t2 = ((-b - math.sqrt(discriminant))/(2*a))
    if t1 > 0 and t2 > 0:
        if t1 < t2:
            t = t1
        elif t1 >= t2:
            t = t2
        return t
    else:
        return None

def ray_from_outside_sphere(normal, d_vec):
    """Checks if a ray comes from outside or inside the sphere.
    Returns True for rays coming from the outside"""
    outside = np.dot(normal, d_vec)
    if outside < 0:
        return True
    if outside > 0:
        return False
    
def reflection(intersection_point, obj, incoming):
    """A function to return a reflected vector for an incoming vector"""
    normal = 1.01*normalize_vector(intersection_point - obj.centerpoint)
    reflected_dir = incoming - 2*(np.dot(incoming, normal))*normal
    return Light(start_point=intersection_point, vec=reflected_dir)

def calculate_color(camera, intersection_point, obj):
    """Calculates the color of a point using the Blinn-Phong model"""
    normal = 1.01*normalize_vector(intersection_point - obj.centerpoint)
    nV = normalize_vector(camera - intersection_point)
    nL = normalize_vector(light_source.point - intersection_point)
    nLnV = np.divide(nV + nL, np.linalg.norm(nV + nL))
    
    i_p_ambient = obj.ambient*light_source.ambient
    i_p_diffuse = obj.diffuse*light_source.diffuse*np.dot(nL, normal) 
    i_p_specular = obj.specular*light_source.specular*((np.dot(normal, nLnV)**(np.divide(obj.shininess, 4))))
    c_p = i_p_ambient + i_p_diffuse + i_p_specular
    for s in range(0, 3):
        if c_p[s] < 0:
            c_p[s] = 0
        if c_p[s] > 1:
            c_p[s] = 1
    return c_p

def rgbtohex(r,g,b):
    """ Turns color in RGB to color in hex code
    https://stackoverflow.com/questions/51591456/can-i-use-rgb-in-tkinter"""
    return f'#{r:02x}{g:02x}{b:02x}'

def move(event):
    """Allows the user to press a point on the sphere to move the light source. Does not move in z direction however."""
    global light_source, light_point
    if (image_width/2) - 40 < event.x < (image_width/2) + 40 and image_height/2 - 40 < event.y < image_height/2 + 40: 
        change_x = (event.x - image_width/2)/image_width
        change_y = (event.y - image_height/2)/image_height
        print(event.x, (event.x - image_width/2), event.y, change_x, change_y, light_point)
        if light_point[0] < 2 or light_point[0] > -2:
            light_point[0] += 10*change_x
        if light_point[1] > plane.centerpoint[1] and (light_point[1] + 1.5*change_y) > (plane.centerpoint[1] + plane.radius):
            light_point[1] -= 10*change_y
        print(light_point)
        light_source = Light_source(point = light_point, ambient=np.array([1, 1, 1]), diffuse = np.array([1, 1, 1]), specular=np.array([1, 1, 1]))
        ray_tracing_algorithm()

class WindowFrame(ttk.Frame):
    """Class for UI"""
    def __init__(self):
        Frame.__init__(self, width=image_width, height=image_height, background="#000000")
        self.master.title("Simple Ray Tracing")
        self.grid(column=0, row=0, sticky=(N, W, E, S))
        self.canvas = Canvas(self, width=image_width, height=image_height, bg="#000000")
        self.img = PhotoImage(width=image_width, height=image_height)
        self.canvas.create_image((image_width/2, image_height/2), image=self.img, state="normal")
        self.canvas.grid(row=1, column=0, sticky=W)
        self.intro_label = Label(self, text="Detta program simulerar ljus \n med hjälp av ray tracing. \n Klicka på sfären för att \n välja position av ljuskällan", foreground="#FFFFFF" , background="#000000", font="TkSmallCaptionFont")
        self.intro_label.grid(row=0, column=0, sticky=(N))
        self.canvas.bind("<Button-1>", move)

class Light:
    """A class to represent physical light rays"""
    def __init__(self, start_point, vec):
        self.start_point = start_point
        self.vec = vec
        
class Sphere:
    """A class to represent a sphere and its properties."""
    def __init__(self, centerpoint, radius, ambient, diffuse, specular, shininess, reflectionco):
        self.centerpoint = centerpoint
        self.radius = radius
        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular
        self.shininess = shininess
        self.reflectionco = reflectionco
    
class Light_source:
    """A class to represent a light_source and its properties."""
    def __init__(self, point, ambient, diffuse, specular):
        self.point = point
        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular

light_point = np.array([0.5, 0.5, 0])
s1 = Sphere(centerpoint=np.array([0, 0.2, -1]), radius=0.4, ambient=np.array([0, 0, 0.6]), diffuse = np.array([0, 0, 0.2]), specular=np.array([1, 1, 1]), shininess=100, reflectionco=0.3)
plane = Sphere(centerpoint=np.array([0, -9000, -1]), radius=9000-0.3, ambient=np.array([0.1, 0.1, 0.1]), diffuse = np.array([0.6, 0.6, 0.6]), specular=np.array([1, 1, 1]), shininess=100, reflectionco=0.1)
light_source = Light_source(point = light_point, ambient=np.array([1, 1, 1]), diffuse = np.array([1, 1, 1]), specular=np.array([1, 1, 1]))
objects = [s1, plane]

def main():
    """Create the GUI and start the main loop."""
    global frame
    frame = WindowFrame()
    frame.after(20, ray_tracing_algorithm())
    mainloop()
    
if __name__ == "__main__":
    main()
