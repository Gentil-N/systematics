
class SampleTime:
       seconds: int
       minutes: int
       hours: int
       days: int
       years: int
       def __init__(self, seconds: int = 0, minutes: int = 0, hours: int = 0, days: int = 0, years: int = 0) :
              self.seconds = seconds
              self.minutes = minutes
              self.hours = hours
              self.days = days
              self.years = years
       def get_total_seconds(self) :
              """
              Get total time in seconds
              """
              seconds_per_minutes = 60
              seconds_per_hours = 60 * seconds_per_minutes
              seconds_per_day = 24 * seconds_per_hours
              seconds_per_year = 365 * seconds_per_day
              return self.years * seconds_per_year + \
                     self.days * seconds_per_day + \
                     self.hours * seconds_per_hours + \
                     self.minutes * seconds_per_minutes + \
                     self.seconds