from django.http import HttpResponse


def index(request):
    return HttpResponse("Starting SKA-Low sensitivity database application")
    