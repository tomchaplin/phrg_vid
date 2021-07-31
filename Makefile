.DEFAULT_GOAL := low

low:
	manim -pql scene.py OmegaDef

high:
	manim -pqh scene.py OmegaDef

play_low:
	mpv media/videos/scene/480p15/OmegaDef.mp4

play_high:
	mpv media/videos/scene/1080p60/OmegaDef.mp4

snap:
	manim -spql scene.py OmegaDef
